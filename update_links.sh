SRC_DIR="a"
OLD_TARGET="/home/user/public_html/dev"
SUB="s/dev/qa/"

find $SRC_DIR -type l \
  -lname "$OLD_TARGET/*" -printf \
  'ln -nsf "$(readlink "%p"|sed $SUB)" "$(echo "%p"|sed $SUB)"\n'\
 > script.sh
