# Generate the case
if [ ! -f ${ID}.absolute.yml ]; then
  # Create an Architecture defined case
  topo.case -name 2H4E2H -architecture 2H.4E.2H
  # Apply corrections to make an absolute case
  topo.absolute -case 2H4E2H.yml -corrections corrections.yml -caseout 2H4E2H.absolute
fi