
Load environment:

    source env.sh


Install (once):

    chmod 755 env.sh install.sh
    ./install.sh


Generate 10000 radiative Bhabha events:

     ./bbbrem/BBBrem 45.6 0.01 10000. 0 0 -t -u -v


Analysis:

     python python/analysis.py