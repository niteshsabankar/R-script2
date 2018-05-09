import subprocess

modules = ['brownA2', 'cyanA2', 'darkolivegreenA2', 'greenyellowA2', 'grey60A2', 'lightyellowA2', 'magentaA2', 'pinkA2', 'purpleA2', 'yellowA2']

for module in modules:
        command = "/w2/xyzhang/src/trinity/Analysis/DifferentialExpression/run_GOseq.pl --genes_single_factor " + module + " --GO_assignments ehux_go.txt --lengths ehux.seq_lens --background background.ids"
        print command
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print "Return Code for", module, "module is", process.returncode