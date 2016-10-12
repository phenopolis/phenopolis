#!/bin/bash

# Backup name
DATE=`date +%Y%m%d`

# Create used directories if they don't exist already
[[ -d /var/lib/phenotips/backups/$DATE ]] || mkdir -p /var/lib/phenotips/backups/$DATE
[[ -d /var/lib/phenotips/next ]] || mkdir -p /var/lib/phenotips/next

# Download the latest release
cd /var/lib/phenotips/next
rm -rf *
wget http://phenotips.org/Download --output-document=phenotips.zip

# Unzip
unzip phenotips.zip
rm phenotips.zip
cd phenotips*

# Determine the version of Tomcat and PhenoTips
TOMCAT=`echo /var/lib/tomcat*`
PHENOTIPS=`echo /var/lib/phenotips/next/phenotips*`

# Copy current configuration
for f in "xwiki.cfg" "xwiki.properties" "web.xml" "struts-config.xml" "hibernate.cfg.xml"
do
  cp /var/lib/tomcat*/webapps/ROOT/WEB-INF/$f webapps/phenotips/WEB-INF/
done
cp /var/lib/tomcat*/webapps/solr/WEB-INF/web.xml webapps/solr/WEB-INF/

# Stop the server
/etc/init.d/tomcat* stop

# Backup the database
#
## Detect where the PhenoTips instance is
if [[ -f ${TOMCAT}/webapps/ROOT/WEB-INF/xwiki.properties ]]
then
  HIBERNATE_CONFIG="${TOMCAT}/webapps/ROOT/WEB-INF/hibernate.cfg.xml"
elif [[ -f ${TOMCAT}/webapps/phenotips/WEB-INF/xwiki.properties ]]
then
  HIBERNATE_CONFIG="${TOMCAT}/webapps/phenotips/WEB-INF/hibernate.cfg.xml"
fi

## Extract the mysql host, DB name, username and password
MHOST=`cat $HIBERNATE_CONFIG | grep 'jdbc:mysql://' | sed -r -e 's/.*\/\/([^\/]+)\/.*/\1/'`
MDB=`cat $HIBERNATE_CONFIG | grep 'jdbc:mysql://' | sed -r -e 's/.*\/\/[^\/]+\/([^?<]+).*/\1/'`
MUSER=`cat $HIBERNATE_CONFIG | grep 'jdbc:mysql://' -A 5 | grep 'connection.username' | sed -r -e 's/.*>([^<]+)<.*/\1/'`
MPASS=`cat $HIBERNATE_CONFIG | grep 'jdbc:mysql://' -A 5 | grep 'connection.password' | sed -r -e 's/.*>([^<]+)<.*/\1/'`

## Compute the backup filename and make sure its parent directory exists
if [[ ! -d ${BACKUPDIR:=/var/lib/phenotips/backups/mysql} ]]
then
  mkdir -p $BACKUPDIR
fi
FILE="${BACKUPDIR}/data-${DATE}.sql"

## Dump the database
mysqldump --events --single-transaction $MDB -u $MUSER -h $MHOST -p$MPASS > $FILE
#
# Done backing up

# Swap directories
mv ${TOMCAT}/webapps/ROOT /var/lib/phenotips/backups/${DATE}/ROOT
mv ${TOMCAT}/webapps/solr /var/lib/phenotips/backups/${DATE}/solr
mv /var/lib/phenotips/solrconfig /var/lib/phenotips/backups/${DATE}/solrconfig
mv ${PHENOTIPS}/webapps/phenotips ${TOMCAT}/webapps/ROOT
mv ${PHENOTIPS}/webapps/solr ${TOMCAT}/webapps/solr
mv ${PHENOTIPS}/solrconfig /var/lib/phenotips/solrconfig

# Fix rights
chown -R tomcat:tomcat /var/lib/phenotips/

# Start the server
/etc/init.d/tomcat* start

echo All done!
