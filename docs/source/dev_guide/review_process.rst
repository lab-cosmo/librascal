.. _review_process:

.. contents::
   :local:

Review process
--------------

This short description should serve as a rough guideline on how we deal with
reviews on pull requests (PR). Some PRs are easy, others are hard. Some require
multiple changes and continous feedback, others can be merged in flawlessly. So
this text should communicate our intent and style and not be a go/no-go
situation.

After successfully adding a new feature, doing maintenance, sorting out a bug or
adding new documentation, in short: contributing to the codebase in any form,
the tests we have employed should all pass. These include tests for
functionality on the C++ as well as the python side but also include checks for
code style as far as they can be automated (see our :ref:`Conding convention
<coding-convention>` for details on the style). We also have a custom target for
make to help with autoformatting C++ as well as Python code.

Further, we provide a PR template on github with several checkmarks to help
start and streamline the process. Please make use of it and get rid of all
unnecessary points, which just globber the PR message. To make it easy for
possible reviewers do not hesitate to contact any of the people invovled in the
development of librascal also beforehand. Some thing are easier to talk about in
person.

So once you have added something to the code base and checked locally if the
checks pass. Checks include the code tests as well as the checks for coding
style as far as they can be checked automatically. You should also have added
tests for your added functionality or adjusted tests, which conflict with a new
design. The next step is to merge or rebase with the master branch. If
necessary, add changes to the code to comply with any other added tests coming
in from the updated master branch if they were not resolved in the merging or
rebase process. Then you are ready to create your PR. Upon creating the PR our
continuous integration spins up and does further checks with different compilers
and compiling modes such as debug.

Now it is time for review. Github provides a button for review requests. Most of
the time the suggested persons for review are a good first address. Add some so that they are notified.

To help speed up the process we follow some rough guidelines so you know what
you can expect from a reviewer's role and we know what we have to do:

 #. If someone starts a review, this person is the go-to person to finish the
    review on the PR and should see it through as fast as possible. This means
    following up the the requested changes and staying part of the
    discussion. And it also includes the approval in the end. As fast as
    possible is often the best way to avoid having to re-merge with the master
    branch if it is changing.

 #. Single comments without actually providing a review are allowed and
    encouraged, but should be used sparsely. They can be really helpful if
    something might be overseen. But if you want changes and involve yourself in
    the process, make a proper review and take the responsibility to see it
    through.

Sometimes there is a lot of commit noise and the master changes and sometimes
the proper order of PR, Review, Requested Changed, remerge, etc. can be
cumbersome. There are no rules to make this easy so please decide on the basis
of the fact that we want a clean compiling master branch. To elaborate on this,
here are a few examples

 #. A PR is approved but there are still changes added afterwards. While this
    can not be avoided completely, it should not be the norm. As the PR
    requester and if the changes are not only of cosmetic nature, as a courtesy
    please ask the approving reviewer again for feedback on the changes.

 #. Changes due to merge or rebase with master because the master branch changed
    can not be avoided. If the PR was approved, master changed and the result
    are conflicts with the PR, this of course has to be solved. Again, as a
    courtesy, please ask for feedback on the changes again because
    merges/rebases can be complicated and we want to avoid a broken master
    branch.
