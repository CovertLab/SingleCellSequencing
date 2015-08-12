# One-time setup

## Python

On Sherlock, run `pyenv local dvv_seq` in your cloned directory

## bash_profile

Add the following to your `$HOME/.bash_profile` (using the appropriate path to your cloned directory):

```bash
export PYTHONPATH="/path/to/SingleCellSequencing:$PYTHONPATH"
```

To import the necessary shared libraries, you will need to execute the following *each time* after you log in to Sherlock (alternatively, you can add it as a line to your `$HOME/.bash_profile`):

```bash
module load dvv_seq
```

## MongoDB

* Create an account at mongolab.com

* Sign in and get to the home screen

* Next to "MongoDB Deployments" you'll see three buttons. Click the one that says "Create new".

* For "Cloud provider", select "amazon web services".  For "Location": "Amazon's US East (Virginia) Region (us-east-1)".

* Under "Plan", choose "Single-node" and select "Sandbox" (...it's free).

* For "Database name" write "single_cell_sequencing".

* Click "Create new MongoDB deployment".

* Back on the home screen, click on "single_cell_sequencing".

* Click on the "Users" tab and then select "Add a database user".

* Choose a username and password.  Note that in fireworks, the password is stored in plaintext.

* Note that the information shown at the top ("To connect using the shell") contains the database hostname and database port (the part after the colon is the port).

## Config files for Fireworks

* Run `python initialize.py`

* Run `lpad reset` and choose `Y`es if prompted


# Usage

## Queue

To queue a set of directories to process, run:

```bash
python fw_queue.py /path/to/directory/containing/libraries
```

To see a list of help options, run:

```bash
python fw_queue.py --help
```

## Run interactively

To run fireworks interactively (after having queued your jobs), run:

```bash
rlaunch rapidfire
```

This is useful for debugging, but you probably don't want to use it in production because it runs serially.

Do not do this on a login node (your tasks are too computationally expensive)

## Run using the SLURM Scheduler

To run using the scheduler:

```bash
qlaunch -r rapidfire --nlaunches infinite --sleep 5
```

This command will run forever until you `Ctrl-C` to kill it once you see that all the output and analysis files have been generated.

`qlaunch` is relatively lightweight, so you can probably get away with running it on a login node.
