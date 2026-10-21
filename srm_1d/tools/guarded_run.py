"""guarded_run -- run a command under a hard memory cap and a wall-clock timeout.

Agents (and people) running long or parallel simulations use this so a runaway
job cannot take over the machine. On 2026-10-13 an unattended verification job
grew a per-step Python list for ~2.4 h until RAM and the page file were
exhausted (20-30+ GB, 99 % of RAM), after the agent that started it had
stopped.

    python -m srm_1d.tools.guarded_run --mem-gb 4 --timeout-min 30 -- python my_probe.py ...

On Windows the command runs inside a Job Object:
  * JobMemoryLimit caps the COMMITTED memory of the whole process tree (the
    command and every process it starts, e.g. multiprocessing workers). An
    allocation past the cap fails, so Python raises MemoryError instead of
    exhausting the machine.
  * KILL_ON_JOB_CLOSE: if this guard process dies for any reason (killed,
    session ended), Windows kills the whole tree. Nothing outlives the guard.
  * The guard kills the tree when the timeout expires.
On POSIX the cap is RLIMIT_AS on the child (inherited by its children) and the
timeout kills the child's process group.

Exit status: the command's own, or 124 on timeout. The last line printed is a
one-line summary with the peak job memory (Windows) and the wall time.
"""
import argparse
import os
import subprocess
import sys
import time

TIMEOUT_EXIT = 124


def _windows_job(mem_bytes):
    """Create a Job Object with a job-wide memory cap and kill-on-close."""
    import ctypes
    from ctypes import wintypes

    k32 = ctypes.WinDLL('kernel32', use_last_error=True)

    class IO_COUNTERS(ctypes.Structure):
        _fields_ = [(n, ctypes.c_ulonglong) for n in (
            'ReadOperationCount', 'WriteOperationCount', 'OtherOperationCount',
            'ReadTransferCount', 'WriteTransferCount', 'OtherTransferCount')]

    class BASIC_LIMIT(ctypes.Structure):
        _fields_ = [('PerProcessUserTimeLimit', ctypes.c_longlong),
                    ('PerJobUserTimeLimit', ctypes.c_longlong),
                    ('LimitFlags', wintypes.DWORD),
                    ('MinimumWorkingSetSize', ctypes.c_size_t),
                    ('MaximumWorkingSetSize', ctypes.c_size_t),
                    ('ActiveProcessLimit', wintypes.DWORD),
                    ('Affinity', ctypes.c_size_t),
                    ('PriorityClass', wintypes.DWORD),
                    ('SchedulingClass', wintypes.DWORD)]

    class EXTENDED_LIMIT(ctypes.Structure):
        _fields_ = [('BasicLimitInformation', BASIC_LIMIT),
                    ('IoInfo', IO_COUNTERS),
                    ('ProcessMemoryLimit', ctypes.c_size_t),
                    ('JobMemoryLimit', ctypes.c_size_t),
                    ('PeakProcessMemoryUsed', ctypes.c_size_t),
                    ('PeakJobMemoryUsed', ctypes.c_size_t)]

    JOB_OBJECT_LIMIT_JOB_MEMORY = 0x00000200
    JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE = 0x00002000
    JobObjectExtendedLimitInformation = 9

    k32.CreateJobObjectW.restype = wintypes.HANDLE
    k32.CreateJobObjectW.argtypes = (ctypes.c_void_p, wintypes.LPCWSTR)
    k32.SetInformationJobObject.argtypes = (
        wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD)
    k32.QueryInformationJobObject.argtypes = (
        wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD,
        ctypes.c_void_p)
    k32.AssignProcessToJobObject.argtypes = (wintypes.HANDLE, wintypes.HANDLE)
    k32.TerminateJobObject.argtypes = (wintypes.HANDLE, wintypes.UINT)

    job = k32.CreateJobObjectW(None, None)
    if not job:
        raise OSError(ctypes.get_last_error(), 'CreateJobObject failed')
    info = EXTENDED_LIMIT()
    info.BasicLimitInformation.LimitFlags = (
        JOB_OBJECT_LIMIT_JOB_MEMORY | JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE)
    info.JobMemoryLimit = int(mem_bytes)
    if not k32.SetInformationJobObject(job, JobObjectExtendedLimitInformation,
                                       ctypes.byref(info), ctypes.sizeof(info)):
        raise OSError(ctypes.get_last_error(), 'SetInformationJobObject failed')

    def assign(proc):
        if not k32.AssignProcessToJobObject(job, int(proc._handle)):
            raise OSError(ctypes.get_last_error(), 'AssignProcessToJobObject failed')

    def kill():
        k32.TerminateJobObject(job, TIMEOUT_EXIT)

    def peak_bytes():
        q = EXTENDED_LIMIT()
        k32.QueryInformationJobObject(job, JobObjectExtendedLimitInformation,
                                      ctypes.byref(q), ctypes.sizeof(q), None)
        return int(q.PeakJobMemoryUsed)

    return assign, kill, peak_bytes


def run(cmd, mem_gb, timeout_min, poll_s=1.0):
    """Run ``cmd`` (list) under the cap; return (exit_code, summary_dict)."""
    mem_bytes = int(mem_gb * 1024 ** 3)
    t0 = time.time()
    if os.name == 'nt':
        assign, kill, peak_bytes = _windows_job(mem_bytes)
        proc = subprocess.Popen(cmd)
        try:
            assign(proc)
        except OSError:
            proc.kill()
            raise
    else:
        import resource
        import signal

        def _limit():
            os.setsid()
            resource.setrlimit(resource.RLIMIT_AS, (mem_bytes, mem_bytes))

        proc = subprocess.Popen(cmd, preexec_fn=_limit)

        def kill():
            os.killpg(proc.pid, signal.SIGKILL)

        def peak_bytes():
            return -1
    timed_out = False
    while True:
        try:
            code = proc.wait(timeout=poll_s)
            break
        except subprocess.TimeoutExpired:
            if timeout_min > 0 and time.time() - t0 > timeout_min * 60.0:
                timed_out = True
                kill()
                proc.wait()
                code = TIMEOUT_EXIT
                break
    summary = dict(exit=code, timed_out=timed_out,
                   wall_s=round(time.time() - t0, 1),
                   peak_job_gb=round(peak_bytes() / 1024 ** 3, 2),
                   mem_cap_gb=mem_gb, timeout_min=timeout_min)
    return code, summary


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Run a command under a memory cap and a wall-clock timeout.')
    ap.add_argument('--mem-gb', type=float, default=4.0,
                    help='cap on the committed memory of the whole process tree (default 4)')
    ap.add_argument('--timeout-min', type=float, default=30.0,
                    help='wall-clock limit in minutes; 0 = none (default 30)')
    ap.add_argument('cmd', nargs=argparse.REMAINDER,
                    help='the command, after --')
    a = ap.parse_args(argv)
    cmd = a.cmd[1:] if a.cmd[:1] == ['--'] else a.cmd
    if not cmd:
        ap.error('no command given (put it after --)')
    code, s = run(cmd, a.mem_gb, a.timeout_min)
    print('[guarded_run] ' + ' '.join(f'{k}={v}' for k, v in s.items()), flush=True)
    return code


if __name__ == '__main__':
    sys.exit(main())
