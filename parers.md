# Run PARERS Interactively

* first get an interactive slurm session using stuart account  
srun --partition=cpu-core-sponsored --account=cpu-stuart-sponsored --ntasks=1 --cpus-per-task=1 --nodes=1 --mem-per-cpu=7500M --time=0-01:00:00 --pty /bin/bash

* open shell inside container  
apptainer shell --bind /data/hps/assoc/ /data/hps/assoc/private/stuart/container/parers.sif

* go to directory with parers.cfg configured for run  
cd <path>

* run pipeline  
python3 /parers/parers

* copy output to lab RSS  
cp -R /data/hps/assoc/private/stuart/data/parers/<output> /data/rss/helens/stuart_k/...



# Run PARERS Non-Interactively

* go to directory with parers.cfg configured for run  
cd <path>

* copy parers.slurm to directory  
cp /data/hps/assoc/private/stuart/parers.slurm ./

* (OPTIONAL) edit parers.slurm as needed to define resources for your job http://gonzo/hpcGuide/SlurmScheduler.html  
nano parers.slurm

* submit the job to the scheduler  
sbatch --mail-user=<YOUR EMAIL> parers.slurm

__Note:__ if you don't want to specify --mail-user you can add it to #SBATCH list at top of your parers.slurm

* you can check on the progress of your job by checking the slurm queue or waiting for emails  
squeue -u <YOUR USER>

__NOTE:__ if nothing shows up in the list, you don't have any active jobs and your job is done.

__WARNING:__ You are now running the parers pipeline non-interactively. This means it will wait until the time has been reached for your job to end for user input.  
* If you think your job is not running correctly, you can cancel the job by using:  
scancel <jobID>

* After your job is done you can check resources consumed if you would like  
sedd <jobID>