#!/bin/bash

## Detect where the PhenoTips instance is
if [[ -f /var/lib/tomcat6/webapps/ROOT/WEB-INF/xwiki.properties ]]
then
  HIBERNATE_CONFIG="/var/lib/tomcat6/webapps/ROOT/WEB-INF/hibernate.cfg.xml"
elif [[ -f /var/lib/tomcat6/webapps/phenotips/WEB-INF/xwiki.properties ]]
then
  HIBERNATE_CONFIG="/var/lib/tomcat6/webapps/phenotips/WEB-INF/hibernate.cfg.xml"
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
if [[ $1 == "undated" ]]
then
  FILE="${BACKUPDIR}/data.sql"
else
  FILE="${BACKUPDIR}/data-$(date +"%F").sql"
fi

## Dump the database
mysqldump --events --single-transaction $MDB -u $MUSER -h $MHOST -p$MPASS > $FILE
