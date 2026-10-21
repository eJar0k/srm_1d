"""srm_1d.tools.guarded_run: the memory cap, the timeout, pass-through exit codes."""
import os
import sys

import pytest

from srm_1d.tools import guarded_run as g

PY = sys.executable
pytestmark = pytest.mark.skipif(os.name != 'nt', reason='Job Object path is Windows-only')


def test_passes_through_exit_code_and_reports_peak():
    code, s = g.run([PY, '-c', 'import sys; sys.exit(3)'], mem_gb=1.0, timeout_min=1)
    assert code == 3 and not s['timed_out']
    assert s['peak_job_gb'] > 0.0


def test_memory_cap_turns_a_runaway_into_memoryerror():
    # A runaway that keeps growing (like a per-step list): under a 0.5 GB cap
    # it must get MemoryError after <= ~0.5 GB, not exhaust the machine. The
    # exit code is the number of 64 MB chunks it obtained.
    src = ('import sys\n'
           'chunks = []\n'
           'try:\n'
           '    while len(chunks) < 100:\n'
           '        chunks.append(bytearray(64 * 1024**2))\n'
           'except MemoryError:\n'
           '    pass\n'
           'sys.exit(len(chunks))\n')
    code, s = g.run([PY, '-c', src], mem_gb=0.5, timeout_min=1)
    assert 1 <= code <= 8          # 8 x 64 MB = 0.5 GB
    assert s['peak_job_gb'] <= 0.55


def test_cap_covers_child_processes():
    # A grandchild is in the same job, so it is capped too.
    inner = ('import sys\n'
             'try:\n'
             '    b = bytearray(2 * 1024**3)\n'
             'except MemoryError:\n'
             '    sys.exit(7)\n')
    src = ('import subprocess, sys\n'
           f'sys.exit(subprocess.call([sys.executable, "-c", {inner!r}]))\n')
    code, _ = g.run([PY, '-c', src], mem_gb=0.5, timeout_min=1)
    assert code == 7


def _process_exited(pid, wait_ms):
    import ctypes
    k32 = ctypes.WinDLL('kernel32', use_last_error=True)
    SYNCHRONIZE = 0x00100000
    h = k32.OpenProcess(SYNCHRONIZE, False, int(pid))
    if not h:
        return True                      # already gone
    try:
        return k32.WaitForSingleObject(h, int(wait_ms)) == 0   # WAIT_OBJECT_0
    finally:
        k32.CloseHandle(h)


def test_killing_the_guard_kills_the_job(tmp_path):
    # Orphan protection: the 2026-10-13 runaway outlived the agent that
    # started it. If the guard dies, Windows must kill its whole job.
    import subprocess
    import time
    pidfile = tmp_path / 'child.pid'
    child = (f'import os, time; open({str(pidfile)!r}, "w").write(str(os.getpid())); '
             'time.sleep(120)')
    guard = subprocess.Popen([PY, '-m', 'srm_1d.tools.guarded_run', '--mem-gb', '1',
                              '--timeout-min', '5', '--', PY, '-c', child])
    t0 = time.time()
    while not pidfile.exists() or not pidfile.read_text():
        assert time.time() - t0 < 30, 'child never started'
        time.sleep(0.1)
    pid = int(pidfile.read_text())
    assert not _process_exited(pid, 0)
    guard.kill()
    guard.wait()
    assert _process_exited(pid, 10000)


def test_timeout_kills_the_job():
    code, s = g.run([PY, '-c', 'import time; time.sleep(60)'], mem_gb=1.0,
                    timeout_min=2.0 / 60.0, poll_s=0.2)
    assert code == g.TIMEOUT_EXIT and s['timed_out']
    assert s['wall_s'] < 30
