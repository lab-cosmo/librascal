.. _review_process:

.. contents::
   :local:

Review process
--------------

This short description should serve as a rough guideline on how we deal with
reviews on pull requests (PR). Some PRs are easy, others are hard. Some require
multiple changes and continuous feedback, others can be merged in
flawlessly. The purpose of this text is to communicate our intent and style and
should not be seen as *law*.

After successfully adding a new feature, doing maintenance, sorting out a bug or
adding new documentation, in short: contributing to the codebase in any form,
you should ensure that all existing and new tests pass. These include tests for
functionality on the C++ as well as the python side but also include tests for
code style as far as they can be automated (see our :ref:`Coding convention
<coding-convention>` for details on the style). We also have a custom ``make``
target to help with autoformatting C++ as well as Python code.

Further, we provide a PR template on Github with several checkmarks to help
start and streamline the review process. Please make use of it, although you can
remove points that are irrelevant to the PR. This template is intended to help
you in preparing your contribution for the review process and includes the main
points most commonly raised by reviewers. Also we highly encourage the use of
the *Draft* PR functionality Github provides. We see it as a way to start a
discussion about the contribution and the proposed changes and it does not have
to be polished at that stage. The difference between a Draft PR and a normal PR
is minor. The Draft stage signals that your contribution is not ready for review
but up for discussion. But it is nevertheless treated by like a normal PR in the
sense that all the continuous integration tools are applied and there is an
immediate feedback from automatic tests. Do not hesitate to contact us at any
stage if you are unsure about anything.

Once you have added something to the code base, please check locally if the
tests pass. These tests include the code tests as well as the tests for code
style as far as they can be checked automatically. You should also add tests for
new functionality. Likewise, if your changes break existing tests, be sure to
update (not deactivate!) the affected tests so that they pass. The next step is
to merge or rebase with the master branch (if you are unsure about which option,
do not hesitate to contact us). If necessary, add changes to the code to comply
with any other added tests coming in from the updated master branch if they were
not resolved in the merging or rebase process. If your PR started out as a
Draft, this is now the time to convert it to a normal PR. If the PR never went
into a Draft PR this is the point where you would actually create the PR. Upon
creating the PR the continuous integration tools run a set of extensive tests
with different compilers and compiling modes, such as debug (mode), to ensure
compatibility with a variety of pre-defined defined sets of options.

Once the PR is ready for review (all items on the checklist are completed), use
the Github button to notify specific developers to review the changes.

To help speed up the process we follow some rough guidelines so everyone knows
what to expect from a reviewer's role:

 * The first person to leave a full review becomes the responsible reviewer for
   that PR.

 * Once the reviewers' concerns have been addressed, the same person is
   responsible for approving them as soon as possible. This not only ensures
   that the changes are available to everyone as soon as they are ready, but
   also minimizes the need for further merges with master and associated extra
   work.

 * Single comments without actually providing a review are allowed and
   encouraged, but should be used sparsely. They can be really helpful if
   something might be overseen. But if you want changes and involve yourself in
   the process, make a proper review and take the responsibility to see it
   through.

 * If you are a reviewer and you start adding changes yourself (excluding the
   one-line propositions that Github provides or just cosmetics) another person
   should become the reviewer. This is to ensure that any part is reviewed by
   someone else other than the creator.

A last note related to our own experience. If it happens that your PR falls into
a busy time with a lot of changes on the master branch which are unrelated to
you, do not hesitate to contact us if the any issues arise *after* your PR has
been approved. While we discourage changes after the PR has been approved,
sometimes it is necessary to be able to merge with master. There are no rules to
make this easy so please decide on the basis of the fact that we want a clean
compiling master branch. To elaborate on this, here are a few examples

 * A PR is approved but there are still changes added afterwards. While this can
   not be avoided completely, it should not be the norm. As the PR requester and
   if the changes are not only of cosmetic nature, as a courtesy please ask the
   approving reviewer again for feedback on the changes.

 * Changes due to merge or rebase with master because the master branch changed
   can not be avoided. If the PR was approved, master changed and the result are
   conflicts with the PR, this of course has to be solved. Again, as a courtesy,
   please ask for feedback on the changes again because merges/rebases can be
   complicated and we want to avoid a broken master branch.
