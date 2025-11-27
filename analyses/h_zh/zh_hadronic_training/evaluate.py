import sys, os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sklearn
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--tag", type=str, default="", help="Input tag")
parser.add_argument("-i", "--input", type=str, default="/ceph/submit/data/group/fcc/ee/analyses/zh/hadronic/training/bdt_model.pkl", help="Input pkl file")
parser.add_argument("-o", "--outDir", type=str, default="/home/submit/jaeyserm/public_html/fccee/h_zh_hadronic/plots_ecm240/training/", help="Output directory")
args = parser.parse_args()


def plot_roc():
    print("Plot ROC")
    train_probs = bdt.predict_proba(train_data)
    train_preds = train_probs[:,1]
    train_fpr, train_tpr, threshold = sklearn.metrics.roc_curve(train_labels, train_preds)
    train_roc_auc = sklearn.metrics.auc(train_fpr, train_tpr)

    test_probs = bdt.predict_proba(test_data)
    test_preds = test_probs[:,1]
    test_fpr, test_tpr, threshold = sklearn.metrics.roc_curve(test_labels, test_preds)
    test_roc_auc = sklearn.metrics.auc(test_fpr, test_tpr)

    # Plot the ROC curve
    plt.figure(figsize=(8, 6))
    plt.plot(train_fpr, train_tpr, color='blue', label=f"Training ROC (AUC = {train_roc_auc:.2f})")
    plt.plot(test_fpr, test_tpr, color='red', label=f"Testing ROC (AUC = {test_roc_auc:.2f})")
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray', label='Random Guess')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend()
    plt.grid()
    plt.savefig(f"{outDir}/roc.png")
    plt.savefig(f"{outDir}/roc.pdf")
    plt.close()


def plot_score():
    print("Plot score")
    train_predictions = bdt.predict_proba(train_data)[:,1]
    test_predictions = bdt.predict_proba(test_data)[:,1]

    # Separate the data into signal and background samples
    train_signal_scores = train_predictions[train_labels == 1]
    train_background_scores = train_predictions[train_labels == 0]
    test_signal_scores = test_predictions[test_labels == 1]
    test_background_scores = test_predictions[test_labels == 0]

    # Plot the BDT scores for signal and background events
    plt.figure(figsize=(8, 6))
    plt.hist(train_signal_scores, bins=50, range=(0, 1), histtype='step', label='Training Signal', color='blue', density=True)
    plt.hist(train_background_scores, bins=50, range=(0, 1), histtype='step', label='Training Background', color='red', density=True)
    plt.hist(test_signal_scores, bins=50, range=(0, 1), histtype='step', label='Testing Signal', color='blue', linestyle='dashed', density=True)
    plt.hist(test_background_scores, bins=50, range=(0, 1), histtype='step', label='Testing Background', color='red', linestyle='dashed', density=True)
    plt.xlabel('BDT Score')
    plt.ylabel('Number of Events (normalized)')
    plt.title('BDT Score Distribution')
    plt.legend()
    plt.grid()
    plt.savefig(f"{outDir}/score.png")
    plt.savefig(f"{outDir}/score.pdf")
    plt.close()

def plot_importance():
    print("Plot importance")
    fig, ax = plt.subplots(figsize=(12, 6))

    importance = bdt.get_booster().get_score(importance_type='weight')
    print(importance)
    sorted_importance = sorted(importance.items(), key=lambda x: x[1], reverse=False)
    print(sorted_importance)
    sorted_indices = [int(x[0][1:]) for x in sorted_importance] # sorted indices
    print(sorted_indices)

    # Get the sorted variable names and their corresponding importances
    sorted_vars = [variables[i] for i in sorted_indices]
    sorted_values = [x[1] for x in sorted_importance]

    # Create a DataFrame and plot the feature importances
    importance_df = pd.DataFrame({'Variable': sorted_vars, 'Importance': sorted_values})
    importance_df.plot(kind='barh', x='Variable', y='Importance', legend=None, ax=ax)
    ax.set_xlabel('BDT score')
    ax.set_title("BDT variable scores", fontsize=16)
    plt.savefig(f"{outDir}/importance.png")
    plt.savefig(f"{outDir}/importance.pdf")
    plt.close()



if __name__ == "__main__":
    baseOutDir = "/home/submit/jaeyserm/public_html/fccee/h_zh_hadronic/plots_ecm240/training/"
    outDir = f"{baseOutDir}/{args.tag}/"
    inputFile = f"/ceph/submit/data/group/fcc/ee/analyses/zh/hadronic/training/{args.tag}.pkl"

    if not os.path.exists(inputFile):
        print(f"Input file {inputFile} does not exist")
        quit()
    
    os.system(f"mkdir -p {outDir}")
    os.system(f"cp {baseOutDir}/index.php {outDir}")

    res = pickle.load(open(inputFile, "rb"))
    bdt = res['model']
    train_data = res['train_data']
    test_data = res['test_data']
    train_labels = res['train_labels']
    test_labels = res['test_labels']
    variables = res['variables']

    plot_score()
    plot_roc()
    plot_importance()