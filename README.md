Setup
=====
Clone this repository into your home directory:

```bash
cd $HOME
git clone https://github.com/CovertLab/SingleCellSequencing.git
```

Note that this code won't work from a share folder in `/home/share`.
This is because `/home/share` isn't mounted on the compute nodes.

Data for a sequenced library should be copied into `$HOME/SingleCellSequencing`

Usage
=====
Suppose you have a `Library_0` sub-directory with the following directory layout:

```bash
$ pwd
$HOME/SingleCellSequencing
$ ls -1 Library_0/
cufflinks_output
tophat_output
```

Then simply run:

```bash
LIBRARY_DIR="Library_0/" make queue-jobs
```

To see how your jobs are doing, run:

```bash
qstat | grep $USER
```

When your jobs are done, you should have files named, for example, `job_cufflinks.sh.o<JOB_ID>` which have the output.
Processed data should appear in the proper folders.

