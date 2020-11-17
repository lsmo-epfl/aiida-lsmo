# aiida-lsmo example submission scripts

This folder contains python scripts for sumitting `aiida-lsmo` workflows.

In order to run the examples, you need
 * a working AiiDA installation
 * the corresponding codes set up in AiiDA (e.g. raspa, zeo++, cp2k, ...)

The examples can be run directly from the command line. For example:
```
./test_multistage_aluminum.py --cp2k-code <cp2k-code-label>
```

 * `run_...` scripts run the workflow in the python interpreter.
   They are designed to run quickly and may not use converged input parameters.
 * `test_...` scripts are like `run_...` scripts but are also part of the continuous integration tests of the plugin (i.e. they are tested on every commit)
 * `submit_...` scripts submit the workflow to the AiiDA daemon and may run for a longer period of time
