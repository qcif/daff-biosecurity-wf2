# Integration tests

These are thin tests, in that they run the entire Python-side workflow without making any assertions - only ensuring that each test case runs without error. The strength of these tests is that each test case covers a different scenario that has raised errors in the past.

Requirements of running this suite:

- Virtual environment at $project_root/venv with requirements installed
- Taxdump must be available to pass P0 validation (specify location in run_tests.sh)
- Bash shell available to run run_tests.sh
