source "firewall-functions.sh"

echo UCL
firewall-setup "128.40.0.0/16"

echo Eduroam
firewall-setup "193.60.0.0/16"

echo Panos home
firewall-setup "89.142.66.56"

echo IoO
firewall-setup "144.82.0.0/16"

echo Sergey
firewall-setup "131.111.80.7"

echo Arcadio
firewall-setup "131.111.186.91"

echo Tom Vulliamy home
firewall-setup "82.16.147.217"

echo Tom Vulliamy work
#firewall-setup "138.37.132.235"
firewall-setup "138.37.0.0/16"

echo  Vincent Plagnol Cambridge
firewall-setup "131.111.84.76"

echo  Nikolas Pontikos home
firewall-setup "79.66.196.129"

echo Mary Fortune Trinity
firewall-remove "131.111.184.8"

echo Arcadio home
firewall-remove "213.205.251.173"

firewall-setup "213.205.251.210"

echo James Poulter
firewall-setup "129.11.65.144"

echo Rahel Gillespie
firewall-setup "130.88.238.38"


