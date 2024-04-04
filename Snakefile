import TreeGenerator
import multiprocessing
import glob
import os

# Dataset
sample_names = [os.path.splitext(os.path.basename(file))[0] for file in glob.glob("dataset/*.fasta")]

# Define rule to generate all desired outputs
rule all:
    input:
        expand("results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.bestTree", sample=sample_names)
        # expand("results/mafft/{sample}/tree/{sample}Tree.svg", sample=sample_names),    
        # expand("results/clustal_omega/{sample}/tree/{sample}Tree.svg", sample=sample_names),
        # expand("results/Muscle/{sample}/tree/{sample}Tree.svg", sample=sample_names),
        # expand("iqtree/{sample}/mafft/mafft.log", sample=sample_names),
        # expand("iqtree/{sample}/clustal_omega/clustal_omega.log", sample=sample_names),
        # expand("iqtree/{sample}/Muscle/Muscle.log", sample=sample_names)

#Alignment rules
rule mafft:
    input:
        data="dataset/{sample}.fasta"
    output:
        "results/mafft/{sample}/{sample}fileAligned.fasta"
    conda:
        "envs/yamlfile.yaml"
    params:
        threads=multiprocessing.cpu_count()
    shell:
        """
        start_time=$(date +%s%3N)
        mafft --auto --thread -{params.threads} {input.data} > {output}
        end_time=$(date +%s%3N)
        execution_time=$((end_time - start_time))
        echo "Mafft = $execution_time milliseconds" > results/mafft/{wildcards.sample}/time.txt
        """

rule Clustal_Omega:
    input:
        previous="results/mafft/{sample}/{sample}fileAligned.fasta",
        data="dataset/{sample}.fasta"
    output:
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta"
    conda:
        "envs/yamlfile.yaml"
    params:
        threads=multiprocessing.cpu_count()
    shell:
        """
        start_time=$(date +%s%3N)
        clustalo -i {input.data} -o {output} --auto --threads={params.threads}
        end_time=$(date +%s%3N)
        execution_time=$((end_time - start_time))
        echo "clustal_omega = $execution_time milliseconds" > results/clustal_omega/{wildcards.sample}/time.txt
        """

rule Muscle:
    input:
        previous="results/clustal_omega/{sample}/{sample}fileAligned.fasta",
        data="dataset/{sample}.fasta"
    output:
        "results/Muscle/{sample}/{sample}fileAligned.fasta"
    conda:
        "envs/yamlfile.yaml"
    params:
        threads=multiprocessing.cpu_count()
    shell:
        """
        start_time=$(date +%s%3N)
        muscle -align {input.data} -output {output} -threads {params.threads}
        end_time=$(date +%s%3N)
        execution_time=$((end_time - start_time))
        echo "Muscle = $execution_time milliseconds" > results/Muscle/{wildcards.sample}/time.txt
        """

# Best AICc model selection rules
rule best_AICc_model_mafft:
    input:
        previous="results/Muscle/{sample}/{sample}fileAligned.fasta",
        data="results/mafft/{sample}/{sample}fileAligned.fasta"
    output:
        "results/mafft/{sample}/{sample}fileAligned.txt",
        "results/mafft/{sample}/{sample}fileAligned.fasta.ckp",
        "results/mafft/{sample}/{sample}fileAligned.fasta.log",
        "results/mafft/{sample}/{sample}fileAligned.fasta.out",
        "results/mafft/{sample}/{sample}fileAligned.fasta.topos",
        "results/mafft/{sample}/{sample}fileAligned.fasta.tree"
    conda:
        "envs/yamlfile.yaml"
    shell:
        "modeltest-ng -i {input.data} -t ml > {output}"  

rule best_AICc_model_clustal_omega:
    input:
        previous="results/mafft/{sample}/{sample}fileAligned.txt",
        data="results/clustal_omega/{sample}/{sample}fileAligned.fasta"
    output:
        "results/clustal_omega/{sample}/{sample}fileAligned.txt",
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.ckp",
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.log",
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.out",
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.topos",
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.tree"
    conda:
        "envs/yamlfile.yaml"
    shell:
        "modeltest-ng -i {input.data} -t ml > {output}"

rule best_AICc_model_Muscle:
    input:
        previous="results/clustal_omega/{sample}/{sample}fileAligned.txt",
        data="results/Muscle/{sample}/{sample}fileAligned.fasta"
    output:
        "results/Muscle/{sample}/{sample}fileAligned.txt",
        "results/Muscle/{sample}/{sample}fileAligned.fasta.ckp",
        "results/Muscle/{sample}/{sample}fileAligned.fasta.log",
        "results/Muscle/{sample}/{sample}fileAligned.fasta.out",
        "results/Muscle/{sample}/{sample}fileAligned.fasta.topos",
        "results/Muscle/{sample}/{sample}fileAligned.fasta.tree"
    conda:
        "envs/yamlfile.yaml"
    shell:
        "modeltest-ng -i {input.data} -t ml > {output}"

# Obtaining best AICc model rules
rule obtain_best_AICc_model_mafft:
    input:
        previous="results/Muscle/{sample}/{sample}fileAligned.txt",
        data="results/mafft/{sample}/{sample}fileAligned.txt"
    output:
        "results/mafft/{sample}/{sample}Model.txt"
    shell:
        "grep 'Best model according to AICc' -A 2 {input.data} | tail -n 1 | sed 's/^.* //' > {output}"

rule obtain_best_AICc_model_clustal_omega:
    input:
        previous="results/mafft/{sample}/{sample}Model.txt",
        data="results/clustal_omega/{sample}/{sample}fileAligned.txt"
    output:
        "results/clustal_omega/{sample}/{sample}Model.txt"
    shell:
        "grep 'Best model according to AICc' -A 2 {input.data} | tail -n 1 | sed 's/^.* //' > {output}"

rule obtain_best_AICc_model_Muscle:
    input:
        previous="results/clustal_omega/{sample}/{sample}Model.txt",
        data="results/Muscle/{sample}/{sample}fileAligned.txt"
    output:
        "results/Muscle/{sample}/{sample}Model.txt"
    shell:
        "grep 'Best model according to AICc' -A 2 {input.data} | tail -n 1 | sed 's/^.* //' > {output}"

# Maximum likelihood tree generation rules
rule maximum_likelihood_tree_step_1_mafft:
    input:
        previous="results/Muscle/{sample}/{sample}Model.txt"
        model="results/mafft/{sample}/{sample}Model.txt",
        msa="results/mafft/{sample}/{sample}fileAligned.fasta"
    output:
        "results/mafft/{sample}/{sample}fileAligned.fasta.raxml.bestTree",
        "results/mafft/{sample}/{sample}fileAligned.fasta.raxml.bestModel",
        "results/mafft/{sample}/{sample}fileAligned.fasta.raxml.mlTrees",
        "results/mafft/{sample}/{sample}fileAligned.fasta.raxml.rba",
        "results/mafft/{sample}/{sample}fileAligned.fasta.raxml.startTree"
    conda:
        "envs/yamlfile.yaml"
    params:
        threads=multiprocessing.cpu_count()
    shell:
        """
        model=$(cat {input.model})
        raxml-ng --msa {input.msa} --model $model --threads {params.threads} --seed 333 --tree pars{{100}},rand{{100}}
        """

rule maximum_likelihood_tree_step_1_clustal_omega:
    input:
        previous="results/mafft/{sample}/{sample}fileAligned.fasta.raxml.bestTree",
        model="results/clustal_omega/{sample}/{sample}Model.txt",
        msa="results/clustal_omega/{sample}/{sample}fileAligned.fasta"
    output:
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.bestTree",
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.bestModel",
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.mlTrees",
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.rba",
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.startTree"
    conda:
        "envs/yamlfile.yaml"
    params:
        threads=multiprocessing.cpu_count()
    shell:
        """
        model=$(cat {input.model})
        raxml-ng --msa {input.msa} --model $model --threads {params.threads} --seed 333 --tree pars{{100}},rand{{100}}
        """

rule maximum_likelihood_tree_step_1_Muscle:
    input:
        previous="results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.bestTree",
        model="results/Muscle/{sample}/{sample}Model.txt",
        msa="results/Muscle/{sample}/{sample}fileAligned.fasta"
    output:
        "results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.bestTree",
        "results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.bestModel",
        "results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.mlTrees",
        "results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.rba",
        "results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.startTree"
    conda:
        "envs/yamlfile.yaml"
    params:
        threads=multiprocessing.cpu_count()
    shell:
        """
        model=$(cat {input.model})
        raxml-ng --msa {input.msa} --model $model --threads {params.threads} --seed 333 --tree pars{{100}},rand{{100}}
        """

# Bootstrap replicates generation rules
rule maximum_likelihood_tree_step_2_mafft:
    input:
        previous="results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.bestTree",
        model="results/mafft/{sample}/{sample}Model.txt",
        msa="results/mafft/{sample}/{sample}fileAligned.fasta"
    output:
        "results/mafft/{sample}/{sample}fileAligned.fasta.raxml.bootstraps",
        "results/mafft/{sample}/{sample}fileAligned.fasta.raxml.log",
        "results/mafft/{sample}/{sample}fileAligned.fasta.raxml.rba"
    conda:
        "envs/yamlfile.yaml"
    params:
        threads=multiprocessing.cpu_count(),
        bootstrap_trees=1000
    shell: 
        """
        model=$(cat {input.model})
        raxml-ng --bootstrap --msa {input.msa} --model $model --threads {params.threads} --seed 333 --bs-trees {params.bootstrap_trees}
        """

rule maximum_likelihood_tree_step_2_clustal_omega:
    input:
        previous="results/mafft/{sample}/{sample}fileAligned.fasta.raxml.bootstraps",
        model="results/clustal_omega/{sample}/{sample}Model.txt",
        msa="results/clustal_omega/{sample}/{sample}fileAligned.fasta"
    output:
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.bootstraps",
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.log",
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.rba"
    conda:
        "envs/yamlfile.yaml"
    params:
        threads=multiprocessing.cpu_count(),
        bootstrap_trees=1000
    shell: 
        """
        model=$(cat {input.model})
        raxml-ng --bootstrap --msa {input.msa} --model $model --threads {params.threads} --seed 333 --bs-trees {params.bootstrap_trees}
        """

rule maximum_likelihood_tree_step_2_Muscle:
    input:
        previous="results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.bootstraps",
        model="results/Muscle/{sample}/{sample}Model.txt",
        msa="results/Muscle/{sample}/{sample}fileAligned.fasta"
    output:
        "results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.bootstraps",
        "results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.log",
        "results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.rba"
    conda:
        "envs/yamlfile.yaml"
    params:
        threads=multiprocessing.cpu_count(),
        bootstrap_trees=1000
    shell: 
        """
        model=$(cat {input.model})
        raxml-ng --bootstrap --msa {input.msa} --model $model --threads {params.threads} --seed 333 --bs-trees {params.bootstrap_trees}
        """
# Support values calculation rules
rule maximum_likelihood_tree_step_3_mafft:
    input:
        previous="results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.bootstraps",
        tree="results/mafft/{sample}/{sample}fileAligned.fasta.raxml.bestTree",
        bootstraps="results/mafft/{sample}/{sample}fileAligned.fasta.raxml.bootstraps"
    output:
        "results/mafft/{sample}/{sample}fileAligned.fasta.raxml.bestTree.raxml.support",
        "results/mafft/{sample}/{sample}fileAligned.fasta.raxml.bestTree.raxml.log"
    conda:
        "envs/yamlfile.yaml"
    params:
        threads=multiprocessing.cpu_count()
    shell:
        "raxml-ng --support --tree {input.tree} --bs-trees {input.bootstraps} --threads {params.threads}"
        
rule maximum_likelihood_tree_step_3_clustal_omega:
    input:
        previous="results/mafft/{sample}/{sample}fileAligned.fasta.raxml.bestTree.raxml.support",
        tree="results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.bestTree",
        bootstraps="results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.bootstraps"
    output:
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.bestTree.raxml.support",
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.bestTree.raxml.log"
    conda:
        "envs/yamlfile.yaml"
    params:
        threads=multiprocessing.cpu_count()
    shell:
        "raxml-ng --support --tree {input.tree} --bs-trees {input.bootstraps} --threads {params.threads}"

rule maximum_likelihood_tree_step_3_Muscle:
    input:
        previous="results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.bestTree.raxml.support",
        tree="results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.bestTree",
        bootstraps="results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.bootstraps"
    output:
        "results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.bestTree.raxml.support",
        "results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.bestTree.raxml.log"
    conda:
        "envs/yamlfile.yaml"
    params:
        threads=multiprocessing.cpu_count()
    shell:
        "raxml-ng --support --tree {input.tree} --bs-trees {input.bootstraps} --threads {params.threads}"

# Tree building rules
rule build_tree_mafft:
    input:
        previous="results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.bestTree.raxml.support",
        tree="results/mafft/{sample}/{sample}fileAligned.fasta.raxml.bestTree.raxml.support"
    output:
        img="results/mafft/{sample}/tree/{sample}Tree.svg"
    run:
        TreeGenerator.TreeGenerator(input.tree, output.img)

rule build_clustal_omega:
    input:
        previous="results/mafft/{sample}/tree/{sample}Tree.svg"
        tree="results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.bestTree.raxml.support"
    output:
        img="results/clustal_omega/{sample}/tree/{sample}Tree.svg"
    run:
        TreeGenerator.TreeGenerator(input.tree, output.img)

rule build_Muscle:
    input:
        previous="results/clustal_omega/{sample}/tree/{sample}Tree.svg"
        tree="results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.bestTree.raxml.support"
    output:
        img="results/Muscle/{sample}/tree/{sample}Tree.svg"
    run:
        TreeGenerator.TreeGenerator(input.tree, output.img)

# Merging all .supports files into one
rule merge_trees:
    input:
        "results/mafft/{sample}/{sample}fileAligned.fasta.raxml.bestTree.raxml.support",
        "results/clustal_omega/{sample}/{sample}fileAligned.fasta.raxml.bestTree.raxml.support",
        "results/Muscle/{sample}/{sample}fileAligned.fasta.raxml.bestTree.raxml.support"
    output:
        output="iqtree/{sample}/{sample}_combined_trees.trees"
    shell:
        "cat {input} > {output}"

# IQTree Test

rule IQTree_Mafft:
    input:
        previous="iqtree/{sample}/{sample}_combined_trees.trees",
        align="results/mafft/{sample}/{sample}fileAligned.fasta",
        tree="iqtree/{sample}/{sample}_combined_trees.trees"
    output:
        ckp="iqtree/{sample}/mafft/mafft.ckp.gz",
        iqtree="iqtree/{sample}/mafft/mafft.iqtree",
        log="iqtree/{sample}/mafft/mafft.log",
        model="iqtree/{sample}/mafft/mafft.model.gz",
        treefile="iqtree/{sample}/mafft/mafft.treefile",
        trees="iqtree/{sample}/mafft/mafft.trees"
    shell:
        """
        iqtree -s {input.align} --trees {input.tree} --test-weight --test-au -n 0 --test 10000 -pre mafft
        mv mafft.ckp.gz {output.ckp}
        mv mafft.iqtree {output.iqtree}
        mv mafft.log {output.log}
        mv mafft.model.gz {output.model}
        mv mafft.treefile {output.treefile}
        mv mafft.trees {output.trees}
        """

rule IQTree_Clustal_Omega:
    input:
        previous="iqtree/{sample}/mafft/mafft.ckp.gz",
        align="results/clustal_omega/{sample}/{sample}fileAligned.fasta",
        tree="iqtree/{sample}/{sample}_combined_trees.trees"
    output:
        ckp="iqtree/{sample}/clustal_omega/clustal_omega.ckp.gz",
        iqtree="iqtree/{sample}/clustal_omega/clustal_omega.iqtree",
        log="iqtree/{sample}/clustal_omega/clustal_omega.log",
        model="iqtree/{sample}/clustal_omega/clustal_omega.model.gz",
        treefile="iqtree/{sample}/clustal_omega/clustal_omega.treefile",
        trees="iqtree/{sample}/clustal_omega/clustal_omega.trees"
    shell:
        """
        iqtree -s {input.align} --trees {input.tree} --test-weight --test-au -n 0 --test 10000 -pre clustal_omega
        mv clustal_omega.ckp.gz {output.ckp}
        mv clustal_omega.iqtree {output.iqtree}
        mv clustal_omega.log {output.log}
        mv clustal_omega.model.gz {output.model}
        mv clustal_omega.treefile {output.treefile}
        mv clustal_omega.trees {output.trees}
        """

rule IQTree_Muscle:
    input:
        previous="iqtree/{sample}/clustal_omega/clustal_omega.ckp.gz",
        align="results/Muscle/{sample}/{sample}fileAligned.fasta",
        tree="iqtree/{sample}/{sample}_combined_trees.trees"
    output:
        ckp="iqtree/{sample}/Muscle/Muscle.ckp.gz",
        iqtree="iqtree/{sample}/Muscle/Muscle.iqtree",
        log="iqtree/{sample}/Muscle/Muscle.log",
        model="iqtree/{sample}/Muscle/Muscle.model.gz",
        treefile="iqtree/{sample}/Muscle/Muscle.treefile",
        trees="iqtree/{sample}/Muscle/Muscle.trees"
    shell:
        """
        iqtree -s {input.align} --trees {input.tree} --test-weight --test-au -n 0 --test 10000 -pre Muscle
        mv Muscle.ckp.gz {output.ckp}
        mv Muscle.iqtree {output.iqtree}
        mv Muscle.log {output.log}
        mv Muscle.model.gz {output.model}
        mv Muscle.treefile {output.treefile}
        mv Muscle.trees {output.trees}
        """