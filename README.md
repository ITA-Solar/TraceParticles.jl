# TraceParticles.jl 
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://ITA-Solar.github.io/TraceParticles.jl/dev
[ci-img]: https://github.com/ITA-Solar/TraceParticles.jl/actions/workflows/test.yml/badge.svg?branch=develop
[ci-url]: https://github.com/ITA-Solar/TraceParticles.jl/actions/workflows/test.yml

| Documentation                    | Continuous integration | 
|    :---:                         |     :----:
|[![][docs-dev-img]][docs-dev-url] | [![][ci-img]][ci-url]  |

TraceParticles.jl provides tools for sampling and simulating charged particles in
an electromagnetic field using
[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).
___
### Git commit message manners
The commit message is mainly for other people, so they should be able to understand it now and six months later. Commit messages cannot be longer than one sentence (line) and should start with a tag identifier (see the end of this section).

Use the imperative form of verbs rather than past tense when referring to changes introduced by the commit in question. For example, "Remove property X", not "Removed property X" and not "I removed...". This tense makes picking, reviewing or reverting commits more readable.

Use following tags for commit messages:

       [DEV] : Code development (including additions and deletions)
       [ADD] : Adding new feature
       [DEL] : Removing files, routines
       [FIX] : Fixes that occur during development, but which have essentially no impact on previous work
       [BUG] : Bug with significant impact on previous work -- `grep`-ing should give restricted list
       [OPT] : Optimisation
       [DBG] : Debugging
       [ORG] : Organisational, no changes to functionality
       [SYN] : Typos and misspellings (including simple syntax error fixes)
       [DOC] : Documentation only
       [REP] : Repository related changes (e.g., changes in the ignore list, remove files)
       [UTL] : Changes in utils
       [CLN] : Clean source code. Remove unneccessary code, change formatting, etc.

Commit message examples:

* "[BUG] Add missing initialisation to tg array"
* "[FIX] Add lowercase casting to params"
* "[CLN] Remove unnecessary allocation for dynamic arrays."
