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

// Establish variables for the matrix
matrix_python = ["3.5", "3.6", "3.7"]

for (python_ver in matrix_python) {

    // (Required) Create a new build configuration
    bc = new BuildConfig()

    // (Required) Give the build configuration a name.
    //            This string becomes the stage header on Jenkins' UI.
    //            Keep it short...
    bc.name = "test-py${python_ver}"
    bc.conda_packages = ["python=${python_ver}"]

    // (Required) Use Linux to execute any stages defined here
    // (Note that Jenkins can only be run with Linux, not MacOSX/Windows)
    bc.nodetype = "linux"

    // (Required) Execute a series of commands to set up the build
    bc.build_cmds = [
        "pip install pytest",
        "python setup.py install",
    ]

    // (Optional) Execute a series of test commands
    bc.test_cmds = [
        "pytest fgscountrate/tests/ --junitxml=results.xml",
        // Add a truly magical command that makes Jenkins work for Python 3.5
        "sed -i 's/file=\"[^\"]*\"//g;s/line=\"[^\"]*\"//g;s/skips=\"[^\"]*\"//g' results.xml",
    ]

    // Add the build to the matrix
    matrix += bc

}

// (Required) Submit the build configurations and execute them in parallel
utils.run(matrix)