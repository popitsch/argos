v=0.0.1-SNAPSHOT
read -p "Build the jars? [y/n]" yn
case $yn in
  [Yy]* ) bin/build-java.sh ${v}; break;;
esac
read -p "Deploy the jars? [y/n]" yn
case $yn in
  [Yy]* ) cp bin/argos-core-${v}*.jar /software/ngm/argos/ ; break;;
esac
read -p "Deploy the code to the cibiv repo? [y/n]" yn
case $yn in
  [Yy]* ) git push ; break;;
esac
