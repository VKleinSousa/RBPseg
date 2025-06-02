import argparse
import textwrap
import logging
import sys
import subprocess
import pandas as pd
import os

# Set up logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def get_data_path():
    """
    Dynamically locate the RBPseg/Data directory.
    This ensures the script works regardless of the package's installation path.

    Returns:
        str: Absolute path to the RBPseg/Data directory.
    """
    try:
        # Find the directory where this script is located
        base_dir = os.path.dirname(os.path.abspath(__file__))

        # Traverse up to the RBPseg base directory
        rbpseg_dir = os.path.dirname(base_dir)

        # Path to the Data directory
        data_path = os.path.join(rbpseg_dir, "Data")

        # Check if the Data directory exists
        if not os.path.exists(data_path):
            raise FileNotFoundError(f"Data directory not found at {data_path}")

        return data_path
    except Exception as e:
        raise RuntimeError(f"Error locating the Data directory: {e}")


def run_foldseek(query_path, database_path, output_path, tmp_dir):
    """
    Run FoldSeek using its CLI.

    Parameters:
        query_path (str): Path to the input query file (e.g., PDB file).
        database_path (str): Path to the FoldSeek database.
        output_path (str): Path to save the output matches.
        tmp_dir (str): Temporary directory for FoldSeek.
    """
    try:
        # Construct the FoldSeek command
        command = [
            "foldseek", "easy-search",
            query_path,
            database_path,
            output_path,
            tmp_dir,
            "--format-output", "query,target,alntmscore,qtmscore,prob,qstart,qend"
        ]

        # Execute the command
        subprocess.run(command, check=True)
        logger.info("FoldSeek run completed successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"FoldSeek failed: {e}")
        sys.exit(1)
    except FileNotFoundError:
        logger.error("FoldSeek is not installed or not in your PATH.")
        sys.exit(1)


def parse_foldseek_results(output_path, mode, region_size=50):
    """
    Parse FoldSeek output and select the top hits for each region of the query (for mode 1).

    Parameters:
        output_path (str): Path to the FoldSeek results file.
        mode (int): 0 for TC classification, 1 for region-specific D-class classification.
        region_size (int): Size of the regions to divide the query into.

    Returns:
        pd.DataFrame: Filtered results or top hits.
    """
    try:
        df = pd.read_csv(output_path, sep="\t", header=None)
        df.columns = ['query', 'target', 'alntmscore', 'qtmscore', 'prob', 'qstart', 'qend']

        filtered_df = df[df['alntmscore'] > 0.4].copy()

        if filtered_df.empty:
            lower_confidence_df = df[df['alntmscore'] > 0.2]
            if lower_confidence_df.empty:
                logger.info("No matches, possibly novel fiber class.")
                return df
            else:
                logger.info("No confident match. But other possible matches (alntmscore > 0.2) are:")
                logger.info(lower_confidence_df)
                return lower_confidence_df

        if mode == 1:
            # Add region column for region-specific grouping
            filtered_df['region'] = ((filtered_df['qstart'] - 1) // region_size) + 1

            # Corrected D_class extraction
            filtered_df['D_class'] = filtered_df['target'].str.extract(r'_D(\d+)_')[0]

            # Group by query and region, then select the top hit for each region
            top_hits = filtered_df.loc[
                filtered_df.groupby(['query', 'region'])['alntmscore'].idxmax()
            ]

            # Sort results for readability
            top_hits = top_hits.sort_values(by=['query', 'region'])

            # Print the top hits for each region
            logger.info("Top hits for each query region:")
            print(top_hits)
            return top_hits

        if mode == 0:
            # Handle TC_class (as in previous implementation)
            filtered_df['TC_class'] = filtered_df['target'].str.extract(r'RBP_\d+_TC_(\d+)')[0]
            consensus_class = (
                filtered_df.groupby('TC_class')['qtmscore']
                .sum()
                .idxmax()
            )
            logger.info(f"Most probable class: {consensus_class}")

            other_classes = filtered_df[filtered_df['TC_class'] != consensus_class]['TC_class'].unique()
            if len(other_classes) > 0:
                logger.info("Possible other classes (lower confidence):")
                logger.info(other_classes)
            else:
                logger.info("No other significant classes.")

        return filtered_df

    except pd.errors.ParserError as e:
        logger.error(f"Error parsing FoldSeek results: {e}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error while parsing FoldSeek results: {e}")
        return None


    except pd.errors.ParserError as e:
        logger.error(f"Error parsing FoldSeek results: {e}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error while parsing FoldSeek results: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
            ___________________________________________
                        RBPseg-classify v0.2.0
            ___________________________________________
            Use this module to classify your input tail fiber model into one of the TC classes or D-classes.
            '''),
        epilog=textwrap.dedent('''\
            _______________________________
            Developed by Victor Klein-Sousa.

            If you used this script, please consider citing us:
            Towards a complete phage tail fiber structure atlas. 
            ______________________________
            ''')
    )

    _defaults = {'target_db': 0}

    parser.add_argument("-p", "--pdb", required=True, help="Path to the input PDB file.")
    parser.add_argument("-o", "--output_path", required=True, help="Directory to store output results.")
    parser.add_argument(
        "-db", "--target_db", type=int, choices=[0, 1], default=_defaults['target_db'],
        help="Classification target: 0 for TC classes, 1 for domain classification."
    )

    args = parser.parse_args()

    logger.info("Starting the main process")
    try:
        os.makedirs(args.output_path, exist_ok=True)
        tmp = os.path.join(args.output_path, 'tmp')
        data_dir = get_data_path()

        # Construct paths to the TC and D databases
        TC_db_path = os.path.join(data_dir, "TC_DB/TC_DB")
        D_db_path = os.path.join(data_dir, "D_DB/D_DB")
        
        

        if args.target_db == 0:
            results_path = os.path.join(args.output_path, "TC_matches.tsv")
            run_foldseek(args.pdb, TC_db_path, results_path, tmp)
            results_df = parse_foldseek_results(results_path, args.target_db)
        elif args.target_db == 1:
            results_path = os.path.join(args.output_path, "D_matches.tsv")
            run_foldseek(args.pdb, D_db_path, results_path, tmp)
            results_df = parse_foldseek_results(results_path, args.target_db)
        else:
            raise ValueError("Invalid target_db mode")

        
        if results_df is not None:
            logger.info("Complete.")
            #logger.info(results_df.head())
        else:
            logger.warning("No results to display.")

    except Exception as e:
        logger.critical(f"An error occurred: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
