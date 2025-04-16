import numpy as np
import torch
import os
import time
import glob

def nsnmf_gpu(V, n_components=10, theta=0.3, max_iter=300, device='cuda'):
    """
    Runs the non-smooth non-negative matrix factorization (nsNMF) algorithm on GPU or CPU.
    
    Parameters:
        V: Input matrix (numpy array)
        n_components: Number of factors (K value)
        theta: Smoothing parameter
        max_iter: Maximum number of iterations
        device: Device to use ('cuda' or 'cpu')
    
    Returns:
        W: Factor matrix W
        H_smoothed: Smoothed factor matrix H (S @ H)
    """
    # Move the input matrix to the specified device
    V = torch.tensor(V, dtype=torch.float32).to(device)
    m, n = V.shape
    r = n_components

    # Randomly initialize W and H (non-negative)
    W = torch.rand((m, r), device=device)
    H = torch.rand((r, n), device=device)

    # Construct the smoothing matrix S (of size r x r)
    I = torch.eye(r, device=device)
    ones = torch.ones((r, r), device=device)
    S = (1 - theta) * I + theta * ones / r

    epsilon = 1e-8  # Prevent division by zero

    # Iteratively update W and H
    for _ in range(max_iter):
        WH = W @ (S @ H)
        W = W * ((V @ H.T @ S.T) / (WH @ H.T @ S.T + epsilon))
        W = torch.clamp(W, min=epsilon)

        WH = W @ (S @ H)
        # Place the smoothing matrix S to the left to ensure the matrix multiplication dimensions match
        H = H * (S @ (W.T @ V) / (S @ (W.T @ WH) + epsilon))
        H = torch.clamp(H, min=epsilon)

    return W.cpu().numpy(), (S @ H).cpu().numpy()

def save_nmf_to_csv(W, H, output_dir='output_csv', prefix='result'):
    """
    Saves the factor matrices W and H to CSV files.
    
    Parameters:
        W, H: Matrices to save
        output_dir: Output directory path
        prefix: File name prefix (suggest to include sample name, K value, and theta info)
    """
    os.makedirs(output_dir, exist_ok=True)
    w_path = os.path.join(output_dir, f'{prefix}_W.csv')
    h_path = os.path.join(output_dir, f'{prefix}_H.csv')

    np.savetxt(w_path, W, delimiter=',')
    np.savetxt(h_path, H, delimiter=',')

    print(f"[✔] Saved W to {w_path}")
    print(f"[✔] Saved H to {h_path}")

if __name__ == "__main__":
    # Check if GPU is available, otherwise use CPU
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Using device: {device}")

    # Use glob to find all scaled_matrix.txt files in the sample directories
    file_list = glob.glob('/mnt/fast00/Minghui/matrix_scale_data/*/scaled_matrix.txt')
    print(f"Found {len(file_list)} scaled_matrix.txt files")

    # Record the runtime for nsNMF on all samples at different K values (format: (sample name, K, runtime))
    runtime_log = []

    # Iterate over each sample file
    for file in file_list:
        # Get the sample name from the directory name containing the file
        sample_name = os.path.basename(os.path.dirname(file))
        print(f"\n===== Processing sample: {sample_name} =====")
        
        # Load the matrix data (assuming data is tab-delimited)
        try:
            V = np.loadtxt(file, delimiter="\t")
        except Exception as e:
            print(f"Failed to load file {file}: {e}")
            continue

        # Run nsNMF for the current sample for K values ranging from 3 to 20
        for k in range(3, 21):
            print(f"\n--- Sample {sample_name}, K={k} ---")
            start_time = time.time()
            
            # Call the nsNMF algorithm (theta and iteration count can be adjusted as needed)
            W, H = nsnmf_gpu(V, n_components=k, theta=0.3, max_iter=2000, device=device)
            
            runtime = time.time() - start_time
            runtime_log.append((sample_name, k, runtime))
            print(f"nsNMF runtime for {sample_name} (K={k}): {runtime:.4f} seconds")
            
            # Save the decomposition result for the current K value to CSV files; file name prefix includes sample name, K, and theta info
            output_prefix = f"{sample_name}_k{k}_theta03"
            save_nmf_to_csv(W, H, output_dir='output_csv', prefix=output_prefix)

    # Write the runtime log for all samples and K values to a log file
    log_path = os.path.join('output_csv', 'batch_runtime_log.csv')
    with open(log_path, 'w') as f:
        f.write("Sample,K,Runtime(s)\n")
        for sample, k, rt in runtime_log:
            f.write(f"{sample},{k},{rt:.4f}\n")
    print(f"\n[✔] Runtime log saved to {log_path}")
