# About this repository's history

This repository was migrated from Bitbucket (`fathomteam/moab`) after Atlassian
retired the Bitbucket issue tracker. Everything below explains how to read what
came across, and what could not come across.

## Read the first line of an issue, not the author avatar

**Every migrated issue, pull request and comment shows the same GitHub account
as its author.** That account is whoever ran the migration. It is not the person
who wrote the text.

GitHub does not allow an API client to set the author of an issue or a comment.
The issue-import API rejects the field outright:

```
POST /repos/{owner}/{repo}/import/issues  ->  422 Unprocessable Entity
"user" is not a permitted key
```

This is deliberate. If it were permitted, any token holder could forge posts
from any account. There is no flag, scope or endpoint that changes it.

So the real author is recorded in the **text** instead. Every migrated item
opens with an attribution line:

```
> Migrated from Bitbucket [issue #102](https://bitbucket.org/...) — *resolved*
> Reported by **[Johannes Probst](https://github.com/jtprobst)** (@jtprobst) on 2019-02-19
```

and every comment opens with one:

```
> **Iulian Grindeanu (@iulian787)** commented on 2023-09-13
```

About 96% of migrated comments name a real GitHub handle. The rest are external
reporters with no GitHub account linked to any commit here, plus one Atlassian
bot. Their names appear as plain text rather than as a guess at a handle,
because a wrong guess would @-mention an unrelated stranger.

**Git history is different and is fully correct.** Commit authorship travels
inside the commit object, so `git log`, `git blame` and the contributor graph
all show the real authors. Only issues and pull requests carry the limitation
above.

## Where things are

| Content | Where it is |
|---|---|
| Bitbucket issues 1-195 | GitHub issues **#1-195**, same numbers |
| Open Bitbucket PRs | real pull requests |
| Declined PRs whose branch survived | real pull requests, closed |
| Merged and other closed PRs | issues titled `[BB PR #N] <title>`, closed |
| Issue attachments and inline images | assets on the `bitbucket-attachments` release |
| Version history | GitHub Releases, one per version tag |

Issue numbers were deliberately kept identical to Bitbucket, so a `#102`
reference written years ago in a comment still points at the right issue.

## Why merged pull requests are issues

A pull request cannot be created in a closed or merged state. It is always
opened from a branch that still exists, and it only becomes merged by actually
being merged. For historical PRs neither condition holds:

* the source branch was deleted on Bitbucket when the PR was merged, so GitHub
  answers `422 {"field": "head", "code": "invalid"}`; or
* the branch survives but is fully merged, so there is no diff and GitHub
  answers `422 "No commits between master and <branch>"`.

Each one is therefore archived as a closed issue carrying the original
description, every comment, the author attribution, the merge-commit SHA
(rewritten to this repository's history) and a link back to Bitbucket.

Declined PRs were never merged, so where the branch survived they *could* be
opened for real and closed. Those are genuine pull requests with diffs.

## What was removed, and why

Large binary test data was stripped from the history so the repository fits
GitHub's limits. Only blobs over 1 MB that are **not** present in the current
`master` or the active refactor branch were removed, so no file you can check
out today is missing. The pack went from roughly 355 MB to 58 MB.

A pull request that depended on a removed file carries a note saying so.

## Commit SHAs changed

Rewriting history changes every commit hash. SHAs quoted inside migrated
descriptions and comments were translated to the new hashes automatically, and
Bitbucket commit links were repointed here. A SHA that could not be resolved
unambiguously was left untouched, so an occasional old hash in old text will not
resolve.

## Not migrated

* **Pull requests opened from forks.** Their branch lives in someone else's
  repository, so GitHub has nothing to open a pull request against. They are
  listed with reconstruction steps in the migration toolkit.
* **Bitbucket's own review state** (approvals, reviewer lists) beyond what
  appears in the comment stream.
