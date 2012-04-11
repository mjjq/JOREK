# .bash_profile

# Get the aliases and functions
if [ -f /etc/bashrc ]; then
    . /etc/bashrc
fi

# User specific environment and startup programs

PATH=$PATH:$HOME/bin

export PATH

export PS1="\[\033[1;34m\][\$(date +%H:%M)][\h:\w]$\[\033[0m\] " 

module purge
module load hdf5-serial mpich2/intel

case `hostname` in
    nashira*)        export ARCH=nashira ;;
    nfs* | manager*) export ARCH=nfs ;;
    login*)          export ARCH=ccamu ;;
    *) 
esac
