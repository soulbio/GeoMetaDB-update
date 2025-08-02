import re
import os
import time
import sqlite3
import tempfile
import argparse
import asyncio
import urllib.request
from Bio import Entrez
from geofetch import Finder
from GEOparse import get_GEO
import datetime
import concurrent.futures

async def fetch_esummary(gse_id):
    """Fetch basic metadata from GEO's esummary endpoint.

    Args:
        gse_id: The GSE ID

    Returns:
        dict: Basic metadata
    """
    data = {
        "title": None,
        "summary": None,
        "overall_design": None,
        "pubmed_id": None,
        "organism": None,
        "type": None,
    }

    try:
        handle = Entrez.esearch(db="gds", term=gse_id)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            return data

        geo_id = record["IdList"][0]
        handle = Entrez.esummary(db="gds", id=geo_id)
        summary = Entrez.read(handle)
        handle.close()

        doc = summary[0]
        data.update(
            {
                "title": doc.get("title"),
                "summary": doc.get("summary"),
                "overall_design": doc.get("overalldesign"),
                "pubmed_id": doc.get("pubmed"),
                "organism": doc.get("taxon"),
                "type": doc.get("gdsType"),
            }
        )
    except Exception as e:
        print(f"Error fetching esummary for {gse_id}: {e}")

    return data

async def download_series_matrix(gse_id, destdir=None, retries=3):
    """Download and extract all series matrix files for a GSE ID.

    Args:
        gse_id: The GSE ID
        destdir: Directory to save the file (defaults to temp dir)
        retries: Number of download attempts

    Returns:
        list: Paths to the downloaded matrix files or empty list if failed
    """
    if destdir is None:
        destdir = tempfile.gettempdir()

    base_url = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:-3]}nnn/{gse_id}/matrix/"
    )

    # Try to get directory listing
    try:
        with urllib.request.urlopen(base_url) as response:
            page = response.read().decode("utf-8")
            # Find all matrix files in the directory
            matrix_files = re.findall(f"{gse_id}.*?_series_matrix.txt.gz", page)
    except:
        # If directory listing fails, try common patterns
        matrix_files = [f"{gse_id}_series_matrix.txt.gz"]

    downloaded_files = []

    for matrix_file in matrix_files:
        url = base_url + matrix_file
        dest_path = os.path.join(destdir, matrix_file)

        for attempt in range(retries):
            try:
                urllib.request.urlretrieve(url, dest_path)
                os.system(f"gunzip -f {dest_path}")
                matrix_path = dest_path[:-3]  # Remove .gz extension
                if os.path.exists(matrix_path):
                    downloaded_files.append(matrix_path)
                    break
            except Exception as e:
                if attempt == retries - 1:
                    print(
                        f"Failed to download matrix {matrix_file} for {gse_id} after {retries} attempts: {e}"
                    )
                continue

    return downloaded_files

def parse_series_matrix(file_paths):
    """Parse GEO series matrix files to extract metadata.

    Args:
        file_paths: List of paths to matrix files or single file path

    Returns:
        dict: Merged metadata from all files with sample-level mapping
    """
    if isinstance(file_paths, str):
        file_paths = [file_paths]

    data = {
        "title": None,
        "summary": None,
        "overall_design": None,
        "type": None,
        "pubmed_id": None,
        "samples": {}  # Dictionary to map sample titles to their metadata
    }

    if not file_paths:
        return data

    for file_path in file_paths:
        if not os.path.exists(file_path):
            continue

        try:
            with open(file_path, "r", encoding="utf-8") as f:
                lines = f.readlines()

            # Collect sample metadata arrays
            sample_geo_accessions = []  # GSM IDs
            sample_titles = []
            sample_source_names = []
            sample_platform_ids = []
            sample_organisms = []
            sample_characteristics = []

            for line in lines:
                if line.startswith("!Series_title") and not data["title"]:
                    data["title"] = (
                        line.split('"')[1]
                        if '"' in line
                        else line.split("!Series_title\t")[1].strip()
                    )
                elif line.startswith("!Series_summary") and not data["summary"]:
                    data["summary"] = (
                        line.split('"')[1]
                        if '"' in line
                        else line.split("!Series_summary\t")[1].strip()
                    )
                elif (
                    line.startswith("!Series_overall_design")
                    and not data["overall_design"]
                ):
                    data["overall_design"] = (
                        line.split('"')[1]
                        if '"' in line
                        else line.split("!Series_overall_design\t")[1].strip()
                    )
                elif line.startswith("!Series_type") and not data["type"]:
                    data["type"] = (
                        line.split('"')[1]
                        if '"' in line
                        else line.split("!Series_type\t")[1].strip()
                    )
                elif line.startswith("!Series_pubmed_id") and not data["pubmed_id"]:
                    data["pubmed_id"] = line.strip().split()[-1]
                    # data["pubmed_id"] = line.strip().split()[1:]
                elif line.startswith("!Sample_geo_accession"):
                    sample_geo_accessions = [
                        x.strip('"') for x in re.findall(r'"([^"]*)"', line)
                    ] or [x.strip() for x in line.split("\t")[1:]]
                elif line.startswith("!Sample_title"):
                    sample_titles = [
                        x.strip('"') for x in re.findall(r'"([^"]*)"', line)
                    ] or [x.strip() for x in line.split("\t")[1:]]
                elif line.startswith("!Sample_source_name_ch1"):
                    sample_source_names = [
                        x.strip('"') for x in re.findall(r'"([^"]*)"', line)
                    ] or [x.strip() for x in line.split("\t")[1:]]
                elif line.startswith("!Sample_platform_id"):
                    sample_platform_ids = [
                        x.strip('"') for x in re.findall(r'"([^"]*)"', line)
                    ] or [x.strip() for x in line.split("\t")[1:]]
                elif line.startswith("!Sample_organism_ch1"):
                    sample_organisms = [
                        x.strip('"') for x in re.findall(r'"([^"]*)"', line)
                    ] or [x.strip() for x in line.split("\t")[1:]]
                elif line.startswith("!Sample_characteristics_ch1"):
                    sample_characteristics = [x.strip('"') for x in re.findall(r'"([^"]*)"', line)] or [x.strip() for x in line.split("\t")[1:]]

            # Map sample metadata using GSM IDs as keys
            max_samples = max(len(sample_titles), len(sample_geo_accessions), 
                            len(sample_source_names), len(sample_platform_ids), len(sample_organisms))
            
            for i in range(max_samples):
                gsm_id = sample_geo_accessions[i] if i < len(sample_geo_accessions) else f"GSM_unknown_{i+1}"
                sample_title = sample_titles[i] if i < len(sample_titles) else f"Sample_{i+1}"
                data["samples"][gsm_id] = {
                    "gsm_id": gsm_id,
                    "title": sample_title,
                    "source_name": sample_source_names[i] if i < len(sample_source_names) else None,
                    "platform_id": sample_platform_ids[i] if i < len(sample_platform_ids) else None,
                    "organism": sample_organisms[i] if i < len(sample_organisms) else None,
                    "characteristics": sample_characteristics[i] if i < len(sample_characteristics) else None
                }
                    
        except Exception as e:
            print(f"Error parsing series matrix file {file_path}: {e}")

    return data

def geoparse_fallback(gse_id):
    """Fallback to GEOparse if matrix and esummary fail.

    Args:
        gse_id: The GSE ID

    Returns:
        dict: Sample metadata
    """
    data = {"samples": {}}

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            gse = get_GEO(geo=gse_id, destdir=tmpdir, silent=True)
            for gsm_name, gsm in list(gse.gsms.items())[:15]:  # Limit to first 15 samples
                data["samples"][gsm_name] = {
                    "gsm_id": gsm_name,
                    "title": gsm.metadata.get("title", [""])[0],
                    "source_name": gsm.metadata.get("source_name_ch1", [""])[0],
                    "platform_id": gsm.metadata.get("platform_id", [""])[0],
                    "organism": gsm.metadata.get("organism_ch1", [""])[0]
                }
    except Exception as e:
        print(f"Error in GEOparse fallback for {gse_id}: {e}")

    return data

async def get_geo_details_using_matrix(gse_id):
    """Get comprehensive GEO details using GEOmetadb first, then matrix files with fallbacks.

    Args:
        gse_id: The GSE ID

    Returns:
        dict: All available metadata
    """
    # Start with empty details
    details = {
        "title": None,
        "summary": None,
        "overall_design": None,
        "pubmed_id": None,
        "samples": {}  # Dictionary to map sample titles to their metadata
    }
    # Step 2: If not in database, try series matrix files
    matrix_paths = await download_series_matrix(gse_id)
    matrix_paths = list(set(matrix_paths))  # Remove duplicates

    if matrix_paths:
        matrix_data = parse_series_matrix(matrix_paths)
        details.update({k: v for k, v in matrix_data.items() if v})

    # Step 2: Supplement with esummary if needed
    if not all(
        details.get(k)
        for k in ["title", "summary", "overall_design", "organism", "type"]
    ):
        esummary_data = await fetch_esummary(gse_id)
        for k, v in esummary_data.items():
            if v and not details.get(k):
                details[k] = v

    # Step 3: Use GEOparse as last resort for missing sample data
    # If no sample metadata found, use GEOparse fallback
    if not details.get("samples") or not any(details["samples"]):
        geoparse_data = geoparse_fallback(gse_id)
        if geoparse_data.get("samples") and not details.get("samples"):
            details["samples"] = geoparse_data["samples"]

    # Cleanup temp files if they exist
    for matrix_path in matrix_paths:
        try:
            if os.path.exists(matrix_path):
                os.remove(matrix_path)
        except:
            pass

    return details


# GSE Metadata Table
def update_gse_table(sqlite_path, gse_id, metadata):
    """Update the gse table with metadata for a specific GSE ID.
    
    Args:
        sqlite_path: Path to the SQLite database
        gse_id: The GSE ID to update
        metadata: Dictionary containing metadata from get_geo_details_using_matrix
    """
    conn = sqlite3.connect(sqlite_path)
    cursor = conn.cursor()
    
    try:
        # Convert pubmed_id to string with ';tab' separator if it's a list
        # pubmed_id = metadata.get('pubmed_id')
        # if isinstance(pubmed_id, list):
        #     pubmed_id = ';tab'.join(str(x) for x in pubmed_id)
        # elif pubmed_id is None:
        #     pubmed_id = None
        # else:
        #     pubmed_id = str(pubmed_id)

        cursor.execute("""
            INSERT INTO gse (gse, title, summary, overall_design, pubmed_id, type)
            VALUES (?, ?, ?, ?, ?, ?)
        """, (
            gse_id,
            metadata.get('title'),
            metadata.get('summary'),
            metadata.get('overall_design'),
            metadata.get('pubmed_id'),
            metadata.get('type')
        ))
        
        conn.commit()
        print(f"Updated gse table for {gse_id}")
        
    except Exception as e:
        print(f"Error updating gse table for {gse_id}: {e}")
        conn.rollback()
    finally:
        conn.close()

# GSM Sample Metadata Table
def update_gsm_table(sqlite_path, gsm_data):
    """Update the gsm table with sample metadata.
    
    Args:
        sqlite_path: Path to the SQLite database
        gsm_data: Dictionary containing GSM sample data
    """
    conn = sqlite3.connect(sqlite_path)
    cursor = conn.cursor()
    
    try:
        for gsm_id, sample_info in gsm_data.items():
            
            # Insert new record
            cursor.execute("""
                INSERT INTO gsm (gsm, title, source_name_ch1, gpl, organism_ch1, characteristics_ch1)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (
                gsm_id,
                sample_info.get('title'),
                sample_info.get('source_name'),
                sample_info.get('platform_id'),
                sample_info.get('organism'),
                sample_info.get('characteristics')
            ))
        
        conn.commit()
        print(f"Updated gsm table for {len(gsm_data)} samples")
        
    except Exception as e:
        print(f"Error updating gsm table: {e}")
        conn.rollback()
    finally:
        conn.close()

# GSE-GSM relationship table
def update_gse_gsm_table(sqlite_path, gse_id, gsm_data):
    """Update the gse_gsm table to link GSE and GSM IDs.
    
    Args:
        sqlite_path: Path to the SQLite database
        gse_id: The GSE ID
        gsm_data: Dictionary containing GSM sample data
    """
    conn = sqlite3.connect(sqlite_path)
    cursor = conn.cursor()
    
    try:
        for gsm_id in gsm_data.keys():
            # Check if the relationship already exists
            cursor.execute("SELECT COUNT(*) FROM gse_gsm WHERE gse = ? AND gsm = ?", (gse_id, gsm_id))
            exists = cursor.fetchone()[0] > 0
            
            if not exists:
                # Insert new relationship
                cursor.execute("""
                    INSERT INTO gse_gsm (gse, gsm)
                    VALUES (?, ?)
                """, (gse_id, gsm_id))
        
        conn.commit()
        print(f"Updated gse_gsm table for {gse_id} with {len(gsm_data)} samples")
        
    except Exception as e:
        print(f"Error updating gse_gsm table for {gse_id}: {e}")
        conn.rollback()
    finally:
        conn.close()

async def process_and_update_gse(sqlite_path, gse_id):
    """Process a single GSE ID and update all relevant tables.
    
    Args:
        sqlite_path: Path to the SQLite database
        gse_id: The GSE ID to process
    """
    try:
        print(f"Processing {gse_id}...")
        
        # Get comprehensive metadata using the existing function
        metadata = await get_geo_details_using_matrix(gse_id)
        
        if not metadata:
            print(f"No metadata found for {gse_id}")
            return
        
        # Update gse table
        update_gse_table(sqlite_path, gse_id, metadata)
        
        # Update gsm table and gse_gsm table if samples exist
        if metadata.get('samples'):
            update_gsm_table(sqlite_path, metadata['samples'])
            update_gse_gsm_table(sqlite_path, gse_id, metadata['samples'])
        
        print(f"Successfully processed {gse_id}")
        
    except Exception as e:
        print(f"Error processing {gse_id}: {e}")

def get_current_geo_details_from_db(sqlite_path):
    conn = sqlite3.connect(sqlite_path)
    cursor = conn.cursor()

    cursor.execute("SELECT gse FROM gse")
    existing_gses = set(row[0] for row in cursor.fetchall())
    conn.close()
    return existing_gses

def _process_single_gse(args):
    """processes a single GSE ID."""
    sqlite_path, gse_id, delay, email = args
    
    Entrez.email = email
    
    try:
        asyncio.run(process_and_update_gse(sqlite_path, gse_id))
        time.sleep(delay)
        return f"Successfully processed {gse_id}"
    except Exception as e:
        return f"Error processing {gse_id}: {e}"

async def update_database_with_missing_gse_ids(sqlite_path, delay=0.5, email="your.email@example.com"):
    start_time = datetime.datetime.now()
    
    # Set Entrez email
    Entrez.email = email
    
    # Get existing GSE IDs from database
    existing_gses = get_current_geo_details_from_db(sqlite_path)
    print(f"Found {len(existing_gses)} existing GSE IDs in database")
    
    # Get all available GSE IDs using geofetch
    gse_obj = Finder()
    all_gse_list = gse_obj.get_gse_all()
    print(f"Found {len(all_gse_list)} total GSE IDs from GEO")
    
    # Find missing GSE IDs
    missing_gse_ids = set(all_gse_list) - existing_gses
    print(f"Found {len(missing_gse_ids)} missing GSE IDs to update")
    
    if not missing_gse_ids:
        print("No missing GSE IDs found. Database is up to date.")
        return

    # Multiprocessing to process missing GSE IDs in parallel
    missing_list = list(missing_gse_ids)[:20]
    print(f"Processing {len(missing_list)} missing GSE IDs using multiprocessing")
    
    # Create argument tuples for each GSE ID
    args_iter = [(sqlite_path, gse_id, delay, email) for gse_id in missing_list]
    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = list(executor.map(_process_single_gse, args_iter))
        
    # Print results
    for result in results:
        print(result)
    
    end_time = datetime.datetime.now()
    total_time = end_time - start_time
    print(f"\nFinished updating database with {len(missing_list)} new GSE IDs")
    print(f"Total time taken: {total_time}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update GEOmetaDB with missing GSE IDs and their metadata.")
    parser.add_argument("--sqlite_path", type=str, required=True, help="Path to the GEOmetadb SQLite database.")
    parser.add_argument("--delay", type=float, default=0.5, help="Delay in seconds between API calls (default: 0.5).")
    parser.add_argument("--email", type=str, default="your.email@example.com", help="Email address for NCBI Entrez API.")
    args = parser.parse_args()
    
    print(f"Using SQLite database at: {args.sqlite_path}")
    print(f"Delay between calls: {args.delay} seconds")
    print(f"Using email: {args.email}")
    
    # Run the async update function
    asyncio.run(update_database_with_missing_gse_ids(
        sqlite_path=args.sqlite_path,
        delay=args.delay,
        email=args.email
    ))