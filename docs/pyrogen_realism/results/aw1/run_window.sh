#!/bin/bash
# AW1: prompt-window comparison (0.1 vs 0.05 s), auto charge, current defaults. Two guarded lanes after a warm-up run.
cd /c/tmp/sv/aw1
export SRM1D_OPENMOTOR_PATH="C:/Users/ejarocki/Documents/Rocketry/Code Stuff/Erosive Burning Solver/openMotor/openMotor"
PY="C:/Users/ejarocki/.pyenv/pyenv-win/versions/3.10.5/python.exe"
OUT=docs/pyrogen_realism/results/aw1
run() { "$PY" -m srm_1d.tools.guarded_run --mem-gb 4 --timeout-min 30 -- "$PY" docs/pyrogen_realism/probes/as1_valid.py "$1" "$OUT/$2" >> "$OUT/$3.log" 2>&1; }
run '{"motor":"chunc","window":0.1}' laneA.jsonl laneA
( for s in '{"motor":"chunc","window":0.05}' '{"motor":"balls","window":0.1}' '{"motor":"balls","window":0.05}'; do run "$s" laneA.jsonl laneA; done ) &
( for s in '{"motor":"haseg","window":0.1}' '{"motor":"haseg","window":0.05}' '{"motor":"pathfin","window":0.1}' '{"motor":"pathfin","window":0.05}'; do run "$s" laneB.jsonl laneB; done ) &
wait
echo done
