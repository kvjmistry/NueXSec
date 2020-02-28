# Script to run project.py --check ana whilst the jobs come back from fermi grid

for i in {1..200}
do
  project.py --xml searchingfornues_run1_dirt_overlay.xml --stage ntuple --checkana
  echo "sleeping for 1 hour"
  sleep 3600
done
