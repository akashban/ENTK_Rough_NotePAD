from radical.entk import Pipeline, Stage, Task, AppManager
import os

# ------------------------------------------------------------------------------
# Set default verbosity
if os.environ.get('RADICAL_ENTK_VERBOSE') == None:
    os.environ['RADICAL_ENTK_REPORT'] = 'True'


hostname = os.environ.get('RMQ_HOSTNAME', 'localhost')
port = os.environ.get('RMQ_PORT', 5672)

if __name__ == '__main__':

    # Create a Pipeline object
    p = Pipeline() 

   # -----------------------------------------------------------------------------------------------------
    # Create a Stage object
    s1 = Stage()
    # This stage will run 5 Dendron-Lipid systems, each having a unique relative concentration
    s1.name = 'Bicelle'

    for cnt in range(5):

        # Create a Task object
        t = Task()
        t.name = 'my-task'        # Assign a name to the task (optional)
        t.executable = '/bin/echo'   # Assign executable to the task
        t.arguments = ['I am task %s in %s'%(cnt, s1.name)]  # Assign arguments for the task executable

        # Add the Task to the Stage
        s1.add_tasks(t)

    # Add Stage to the Pipeline
    p.add_stages(s1)

   # -----------------------------------------------------------------------------------------------------

    # Need to add an analysis stage here. The objective of this stage would be to distinguish between a pancake like bicelle and 
    # a sinosoidal bicelle.
    # Lets assume 3 of the initial 5 structures are determined to be stable

    # Create another Stage object
    s2 = Stage()
    # This stage will continue running the 3 stable strctures 
    
    s2.name = 'Vesicle'

    for cnt in range(3):

        # Create a Task object
        t = Task()
        t.name = 'my-task'        # Assign a name to the task (optional, do not use ',' or '_')
        t.executable = '/bin/echo'   # Assign executable to the task
        t.arguments = ['I am task %s in %s'%(cnt, s2.name)]  # Assign arguments for the task executable

        # Add the Task to the Stage
        s2.add_tasks(t)

    # Add Stage to the Pipeline
    p.add_stages(s2)

   # -----------------------------------------------------------------------------------------------------

   # Need to add an analysis stage here. The objective of this stage would be to distinguish between a regular and an irregular vesicle.
   # Lets assume only 1 assembly to be a regular vesicle.
  
 # Create another Stage object
    s3 = Stage()
    s3.name = 'Is vesicle Stable'
    # This stage will exceute the production run for 1 Dendron-Lipid assembly. This assembly has been given the go-ahead sign by
    # both the analysis stages.

    for cnt in range(1):

        # Create a Task object
        t = Task()
        t.name = 'my-task'        # Assign a name to the task (optional, do not use ',' or '_')
        t.executable = '/bin/echo'   # Assign executable to the task
        t.arguments = ['I am task %s in %s'%(cnt, s3.name)]  # Assign arguments for the task executable

        # Add the Task to the Stage
        s3.add_tasks(t)

    # Add Stage to the Pipeline
    p.add_stages(s3)

 # -----------------------------------------------------------------------------------------------------

    # Create Application Manager
    appman = AppManager(hostname=hostname, port=port)

    # Create a dictionary describe four mandatory keys:
    # resource, walltime, and cpus
    # resource is 'local.localhost' to execute locally
    res_dict = {

        'resource': 'local.localhost',
        'walltime': 10,
        'cpus': 1
    }

    # Assign resource request description to the Application Manager
    appman.resource_desc = res_dict

    # Assign the workflow as a set or list of Pipelines to the Application Manager
    # Note: The list order is not guaranteed to be preserved
    appman.workflow = set([p])

    # Run the Application Manager
    appman.run()
