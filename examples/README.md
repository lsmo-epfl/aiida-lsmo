# aiida-lsmo example submission scripts

This folder contains python scripts for sumitting examples of `aiida-lsmo` workflows.

For example
```
./test_multistage_aluminum.py --cp2k-code <cp2k-code-label>
```

 * `run_...` scripts run the workflow in the python interpreter.
   They are designed to run quickly (and may not use converged input parameters)
 * `test_...` scripts are like `run_...` scripts but are also as part of the continuous integration tests (i.e. they are tested on every commit)
 * `submit_...` scripts submit the workflow to the AiiDA daemon (e.g. because they run for a longer period of time)
