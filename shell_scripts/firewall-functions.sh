# sudo bash

# yum install httpd

# service httpd start
# chkconfig httpd oo
# chkconfig httpd on 


# remove ports
function remove-ports() {
	firewall-cmd --zone=public --remove-port=80/tcp --permanent
	firewall-cmd --zone=public --remove-port=443/tcp --permanent
	firewall-cmd --zone=public --remove-port=5000/tcp --permanent
	firewall-cmd --zone=public --remove-port=8000/tcp --permanent
	firewall-cmd --zone=public --remove-port=8888/tcp --permanent
}

function add-ports() {
	firewall-cmd --zone=public --add-port=80/tcp
	firewall-cmd --zone=public --add-port=80/tcp --permanent 
	firewall-cmd --zone=public --add-port=443/tcp
	firewall-cmd --zone=public --add-port=443/tcp --permanent 
	firewall-cmd --zone=public --add-port=8080/tcp
	firewall-cmd --zone=public --add-port=8080/tcp --permanent 
	firewall-cmd --zone=public --add-port=8000/tcp
	firewall-cmd --zone=public --add-port=8000/tcp --permanent 
	firewall-cmd --zone=public --add-port=5000/tcp
	firewall-cmd --zone=public --add-port=5000/tcp --permanent
	firewall-cmd --zone=public --add-port=8888/tcp
	firewall-cmd --zone=public --add-port=8888/tcp --permanent
	firewall-cmd --zone=public --add-port=9042/tcp
	firewall-cmd --zone=public --add-port=9042/tcp --permanent
}

function firewall-setup() {
	firewall-cmd --zone=public  --permanent --add-rich-rule="rule family="ipv4" source address="$1" port protocol="tcp" port="80" accept"
	firewall-cmd --zone=public  --permanent --add-rich-rule="rule family="ipv4" source address="$1" port protocol="tcp" port="443" accept"
	firewall-cmd --zone=public  --permanent --add-rich-rule="rule family="ipv4" source address="$1" port protocol="tcp" port="8080" accept"
	firewall-cmd --zone=public  --permanent --add-rich-rule="rule family="ipv4" source address="$1" port protocol="tcp" port="8000" accept"
	firewall-cmd --zone=public  --permanent --add-rich-rule="rule family="ipv4" source address="$1" port protocol="tcp" port="5000" accept"
	firewall-cmd --zone=public  --permanent --add-rich-rule="rule family="ipv4" source address="$1" port protocol="tcp" port="8888" accept"
	firewall-cmd --reload
	firewall-cmd --list-all-zones
	#firewall-cmd --status
	firewall-cmd --get-active-zones
}

function firewall-remove() {
	firewall-cmd --zone=public  --permanent --remove-rich-rule="rule family="ipv4" source address="$1" port protocol="tcp" port="80" accept"
	firewall-cmd --zone=public  --permanent --remove-rich-rule="rule family="ipv4" source address="$1" port protocol="tcp" port="443" accept"
	firewall-cmd --zone=public  --permanent --remove-rich-rule="rule family="ipv4" source address="$1" port protocol="tcp" port="8080" accept"
	firewall-cmd --zone=public  --permanent --remove-rich-rule="rule family="ipv4" source address="$1" port protocol="tcp" port="8000" accept"
	firewall-cmd --zone=public  --permanent --remove-rich-rule="rule family="ipv4" source address="$1" port protocol="tcp" port="5000" accept"
	firewall-cmd --zone=public  --permanent --remove-rich-rule="rule family="ipv4" source address="$1" port protocol="tcp" port="8888" accept"
	firewall-cmd --reload
	firewall-cmd --list-all-zones
	#firewall-cmd --status
	firewall-cmd --get-active-zones
}




