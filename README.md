# GEOmetaDB Database Update Tool

A Python utility developed by **SoulBio** to update the GEOmetaDB SQLite database with missing GEO Series (GSE) records, their associated samples (GSM), and relevant metadata.

This project aims to help researchers work with more up-to-date GEO metadata while maintaining compatibility with the existing GEOmetaDB database structure.

## Download Updated GEOmetaDB

For researchers who need an updated database without running the update script, we've made our latest SQLite database available for download.

**[Download the updated GEOmetaDB SQLite database (Google Drive)](https://drive.google.com/file/d/13rk1EiEItFYqtLN8mEOb1KudL2fSDLlQ/view?usp=sharing)**

You can use this database directly with SQLite-compatible tools and existing GEOmetaDB workflows.

Alternatively, you can use the update script below to update an existing database yourself.

## Features

The tool updates three main tables in the GEOmetaDB SQLite database.

### 1. `gse` — GEO Series Metadata

Adds missing GSE records with the following metadata:

- `title`: Series title
- `summary`: Series description
- `overall_design`: Overall experimental design
- `pubmed_id`: Associated PubMed ID
- `type`: Series type

### 2. `gsm` — GEO Sample Metadata

Adds metadata for associated GSM samples:

- `title`: Sample title
- `source_name_ch1`: Sample source
- `gpl`: Associated platform ID
- `organism_ch1`: Organism information
- `characteristics_ch1`: Sample characteristics

### 3. `gse_gsm` — Series–Sample Relationships

Maintains relationships between GEO Series and their associated samples:

- `gse`: GSE Series ID
- `gsm`: GSM Sample ID

## Usage

### 1. Prepare the SQLite Database

Download an existing GEOmetaDB SQLite database or use the updated database provided above.

We recommend keeping a backup before running the update script.

### 2. Run the Update Script

```bash
python GEOmetaDB_update.py --sqlite_path /path/to/GeoMetaDB.sqlite
```

### Optional Parameters

| Parameter | Description | Default |
|---|---|---|
| `--sqlite_path` | Path to the GEOmetaDB SQLite database | Required |
| `--batch_size` | Number of GSE IDs processed per batch | `10` |
| `--delay` | Delay in seconds between API calls to reduce rate limiting | `1.0` |

### Examples

**Basic usage:**

```bash
python GEOmetaDB_update.py --sqlite_path ./GeoMetaDB.sqlite
```

**Custom batch size and API delay:**

```bash
python GEOmetaDB_update.py \
  --sqlite_path ./GeoMetaDB.sqlite \
  --batch_size 5 \
  --delay 2.0
```

## Contributing

We welcome feedback, bug reports, and contributions from the research and open-source communities.

If you encounter issues with metadata completeness, database compatibility, or the update process, please [open an issue](https://github.com/soulbio/GeoMetaDB-update/issues) or submit a pull request.

## Maintained by SoulBio

Developed by [SoulBio](https://github.com/soulbio) to support open-source bioinformatics tools and improve access to public biological research data.
