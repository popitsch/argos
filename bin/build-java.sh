#!/bin/bash
if [ $# -lt 1 ]; then
  echo "usage: build.sh <version>"
  echo the version should match the entry in the pom.xml file.
  head pom.xml
  exit 1
fi
export MAVEN_OPTS=-Xss2m
mvn $2 clean assembly:assembly -DskipTests=true -U -Denv=dev -Dmaven.wagon.http.ssl.insecure=true -Dmaven.wagon.http.ssl.allowall=true
cp target/argos-core-$1-jar-with-dependencies.jar bin/argos-core-$1.jar
cp target/argos-core-$1-tests.jar bin/argos-core-$1-tests.jar
cp target/argos-core-$1.jar bin/argos-core-$1-library.jar
