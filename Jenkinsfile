// An example Jenkinsfile

// scm_checkout() does the following:
//     1. Disables pipeline execution if [skip ci] or [ci skip] is present in the
//          commit message, letting users exclude individual commits from CI
//     2. Clones the Git repository
//     3. Creates a local cache of the repository to avoid commit drift between tasks
//        (i.e. Each stage is guaranteed to receive the same source code regardless of
//          commits taking place after the pipeline has started.)
if (utils.scm_checkout()) return

// Set up the matrix of builds
matrix = []

// (Required) Create a new build configuration
bc0 = new BuildConfig()

// (Required) Give the build configuration a name.
//            This string becomes the stage header on Jenkins' UI.
//            Keep it short...
bc0.name = "test"
bc.conda_packages = ["python=3.6"]

// (Required) Use Linux to execute any stages defined here
// (Note that Jenkins can only be run with Linux, not MacOSX/Windows)
bc0.nodetype = "linux"

// (Required) Execute a series of commands to set up the build
bc0.build_cmds = [
    "pip install pytest",
    "python setup.py install",
]

// (Optional) Execute a series of test commands
bc0.test_cmds = ["pytest fgscountrate/tests/ --junitxml=results.xml"]

// Add the build to the matrix
matrix += bc0

// (Required) Submit the build configurations and execute them in parallel
utils.run(matrix)