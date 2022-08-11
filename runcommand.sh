#!/bin/bash
#SBATCH --job-name="FULL_TOPOLY_SIKORA"
#SBATCH --gres=gpu
source $HOME/topology_AI_pipeline/venv/bin/activate
python $HOME/topology_AI_pipeline/pipeline.py -i $HOME/topology_AI_pipeline/test_full.txt