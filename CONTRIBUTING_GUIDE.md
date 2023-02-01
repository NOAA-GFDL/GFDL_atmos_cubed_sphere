# Contributing Guide for the GFDL_atmos_cubed_sphere Repository

This guide will walk a developer through the correct process for making a Pull Request in this repository. 
GFDL_atmos_cubed_sphere is a repository with 4 different development branches.  It is imperitive that the developer understands the difference between each branch.

## Understanding the Development Branches

There are 4 different development branches being supported in this repository.

| Branch | Description |
| :--- | :--- |
| main | This branch is the main development branch. The SHiELD model will compile with this branch. When there is a Public Release of the FV3 Dyncamical Core, updates will first be introduced to this branch. |
| dev/gfdl | This branch is used for all AM4 based GFDL Models. |
| dev/emc | This branch is used for the UFS Weather Model development. |
| dev/gfdl_am5 | This branch is being used for GFDL AM5 development. |

## How to contribute code changes

1. Create a Fork
    1. Click on **Fork** in the top right of the repository GitHub page
    2. The **Owner** should be set to your GitHub username
    3. The **Repository Name** should be GFDL_atmos_cubed_sphere
    4. Click **Create fork**

2. Create an Issue describing the change that you would like to implement.
    1. Navigate to the **Issue** tab at the top of the repository GitHub page
    2. Click on the **New issue** button
    3. Choose from one of the suggested templates (Bug Report, Feature Request, or Support Request)
    4. Fill out the Issue with specifics and submit issue

3. Clone the repository locally on your machine

    `git clone https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere.git`

4. Add your fork locally

    This guide will refer to the fork as `myFork`, but you can name this anything.

    `git remote add myFork https://github.com/<username>/GFDL_atmos_cubed_sphere.git`
    
    `git remote -v` will display all remote repositories that you have added.  The repository that you cloned will be named `origin` by default.  
    
    The ouput of `git remote -v` should be similar to:
    
    ```
    myFork https://github.com/<username>/GFDL_atmos_cubed_sphere.git (fetch)
    myFork https://github.com/<username>/GFDL_atmos_cubed_sphere.git (push)
    origin https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere.git (fetch)
    origin https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere.git (push)
    ```

5.  Checkout the branch that you would like your changes added to

    Refer to section [Understanding the Development Branches](#Understanding-the-Development-Branches) to choose which branch to checkout
    
    This guide will reference this branch as `baseBranch`

    `git checkout baseBranch `

6. Create a feature branch to make your changes to

    This guide will refer to this new branch as `newBranch`, but you should name this branch with a unique name related to the task

    `git checkout -b newBranch`

7. Update the code
    1. To see the files you have modified use the command `git status`
    2. To see exactly what you changed in each file use the command `git diff`
    3. When you are satisfied with your changes stage them for a commit

        `git add <filenames or paths of all added and modified files>`

    4. Make a commit
    
        `git commit -m "Descriptive message describing what you have changed in this commit"`
        
    5. Make sure branch is up to date with the base branch (main, dev/gfdl, dev/emc, or dev/gfdl_am5)

        `git fetch origin baseBranch`
        
        `git merge origin/baseBranch`

    6. Push that commit to your fork

        `git push myFork newBranch`
        
8. Create a Pull Request
    1. Navigate to your fork on GitHub

       The URL to get you to your fork should be `https://github.com/<username>/GFDL_atmos_cubed_sphere`
       
    2. Navigate to the **Pull requests** tab at the top of the repository GitHub page
    3. Click on the **New pull request** button
    4. The **base repository** should be *NOAA-GFDL/GFDL_atmos_cubed_sphere*
    5. The **base** branch is the branch you would like to add your changes to

        Refer to section [Understanding the Development Branches](#Understanding-the-Development-Branches)
        
        This is the same branch that you originally checked out in Step 5 of this guide that was referred to as `baseBranch`
    
    6. The **head repository** should be your fork (e.g. *\<username\>/GFDL_atmos_cubed_sphere*)
    7. The **compare** branch is the feeature branch containing your updates.  This was referred to as `newBranch` in this guide
 
        You should now see a comparison of the two branches  
  
    8. Click on the **Create pull request** button
    9. Fill in the details of the Pull request, being sure to follow the template provided:
        1. Provide a desciption: Include a summary of the change and which issue is fixed. 
        Please also include relevant motivation and context. 
        List any dependencies that are required for this change.
        2. Link the Issue from Step 2 by including a line `Fixes #123` where *123* is the Issue #
        3. Please describe the tests that you ran to verify your changes. 
        Please also note any relevant details for your test configuration (e.g. compiler, OS). 
        Include enough information so someone can reproduce your tests.
        4. Ensure that all checkboxes  are populated.  If something does not apply, please check it and note that it does not apply somewhere in the PR.
        If you have not completed an item in the checklist, the PR will not be merged.
        
            To check a box replace the space in `[ ]` with an x `[x]`
        5. Click on the **Create pull request** button.
    10.  Code managers will assign reviewers to the PR. 
    If you would like someone specific to review your PR please leave a comment on the PR requesting that.
    When all reviewers approve the code, a code manager will merge the code and your changes will now be in the relevant development branch.
