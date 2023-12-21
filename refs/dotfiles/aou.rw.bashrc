# .bash_profile for AoU RW
# Updated December 2023
# Ryan Collins <Ryan_Collins@dfci.harvard.edu>

# Colorize shell and utilities
export CLICOLOR=1
export LS_COLORS='ExFxCxDxBxegedabagacad'

# Simple ls with options
alias l='ls -lhtr'

# Default to zless instead of less
alias 'less=zless'

# Colorize grep output
alias 'grep=grep --color=auto'
alias 'egrep=egrep --color=auto'
alias 'fgrep=fgrep --color=auto'

# Always default to using cromshell-alpha
alias 'cromshell=cromshell-alpha'
