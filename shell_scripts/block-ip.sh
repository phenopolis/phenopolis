
function block() {
    iptables -A INPUT -s $1 -j DROP
}

# CSS Certificate Spider (http://www.css-security.com/certificatespider/)
block 137.116.71.170
# shodan
block 198.20.87.98
# somebody in France?
block 37.187.114.171
# scan-15.shadowserver.org
block 184.105.247.196

#service iptables save
/sbin/iptables-save
# /usr/libexec/iptables/iptables.init save

