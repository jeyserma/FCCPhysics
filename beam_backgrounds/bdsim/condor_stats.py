import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

mpl.use("Agg")

def main():

    clusterId = sys.argv[1]
    filename = "/tmp/condor_stats.txt"
    outDir = "/home/submit/jaeyserm/public_html/condor"
    cmd = f"condor_history -constraint ClusterId=={clusterId} -af:V ClusterId ProcId Owner RequestMemory MemoryUsage RemoteWallClockTime RemoteUserCpu RemoteSysCpu RequestCpus ExitCode ExitBySignal ExitSignal > {filename}"
    os.system(cmd)
    

    data_mem = []
    data_cpu = []
    data_time = []

    with open(filename) as f:
        for line in f:
            
            parts = line.split()
            if len(parts) < 6:
                continue

            exitcode = parts[9]
            if exitcode.isdigit() and  int(exitcode) != 0:
                continue

            print(line)
            mem_req = float(parts[3])
            mem_usage = float(parts[4])
            data_mem.append(mem_usage)

            time_wallclock = float(parts[5])
            time_usercpu = float(parts[6])
            time_syscpu = float(parts[7])
            ncores = float(parts[8])

            cpu_eff = (time_usercpu+time_syscpu) / (time_wallclock*ncores)  * 100.
            data_cpu.append(cpu_eff)
            data_time.append(time_wallclock / 3600)


    # plot memory
    NBINS = 60
    plt.figure(figsize=(7,5))

    vmin, vmax = min(data_mem), max(data_mem)
    bins = np.linspace(vmin, vmax, NBINS + 1)

    mean = sum(data_mem) / len(data_mem)
    plt.hist(data_mem, bins=bins, histtype='step', linewidth=1, color='black')

    plt.xlabel("Memory efficiency (mb)")
    plt.ylabel("Counts")
    plt.title(f"Cluster ID {clusterId}, average = {mean:.1f} mb ")

    plt.tight_layout()
    plt.savefig(f"{outDir}/memory.png", dpi=200)
    plt.close()


    # cpu eff
    NBINS = 60
    plt.figure(figsize=(7,5))

    vmin, vmax = min(data_cpu), max(data_cpu)
    bins = np.linspace(vmin, vmax, NBINS + 1)

    mean = sum(data_cpu) / len(data_cpu)
    plt.hist(data_cpu, bins=bins, histtype='step', linewidth=1, color='black')

    plt.xlabel("CPU efficiency (%)")
    plt.ylabel("Counts")
    plt.title(f"Cluster ID {clusterId}, average = {mean:.1f} % ")

    plt.tight_layout()
    plt.savefig(f"{outDir}/cpu.png", dpi=200)
    plt.close()


    # plot time
    NBINS = 60
    plt.figure(figsize=(7,5))

    vmin, vmax = min(data_time), max(data_time)
    bins = np.linspace(vmin, vmax, NBINS + 1)

    mean = sum(data_time) / len(data_time)
    plt.hist(data_time, bins=bins, histtype='step', linewidth=1, color='black')

    plt.xlabel("Job time (h)")
    plt.ylabel("Counts")
    plt.title(f"Cluster ID {clusterId}, average = {mean:.1f} hours ")

    plt.tight_layout()
    plt.savefig(f"{outDir}/time.png", dpi=200)
    plt.close()



if __name__ == "__main__":
    main()