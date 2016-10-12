#UCL
alias life='ssh -Y rmhanpo@life.biol.ucl.ac.uk'
alias research-data='ssh -Y rmhanpo@ssh.rd.ucl.ac.uk'
alias sftp-research-data='sftp rmhanpo@ssh.rd.ucl.ac.uk'
alias pchuckle='ssh -Y rmhanpo@pchuckle.cs.ucl.ac.uk'
alias elwood='ssh rmhanpo@elwood.cs.ucl.ac.uk'
alias mount-pchuckle='ssh -f rmhanpo@life.biol.ucl.ac.uk -L2222:pchuckle.cs.ucl.ac.uk:22 -N ; sshfs -p 2222 rmhanpo@localhost:/home/rmhanpo /Users/pontikos/ucl'
alias mount-IBDAJE='ssh -f rmhanpo@life.biol.ucl.ac.uk -L2223:pchuckle.cs.ucl.ac.uk:22 -N ; sshfs -p 2223 rmhanpo@localhost:/cluster/project8/IBDAJE /Users/pontikos/IBDAJE'
alias zeppo-1='ssh zeppo-10-19'
alias zeppo-2='ssh zeppo-10-20'

function cannon() {
ssh -t cannon "cd $PWD; exec \$SHELL --login"
}


alias phenotips='ssh phenotips.cs.ucl.ac.uk'

