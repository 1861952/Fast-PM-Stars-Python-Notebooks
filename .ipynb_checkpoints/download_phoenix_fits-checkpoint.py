import os
import requests
from bs4 import BeautifulSoup

# URL of the PHOENIX FITS files
base_url = "https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/phoenix/phoenixm40/"

# Where to save the files (customize this if you want)
download_folder = os.path.join(os.getcwd(), "All_Phoenix_Models")
os.makedirs(download_folder, exist_ok=True)

# Fetch the page
response = requests.get(base_url)
soup = BeautifulSoup(response.text, 'html.parser')

# Loop through links and download all .fits files
for link in soup.find_all('a'):
    href = link.get('href')
    if href and href.endswith('.fits'):
        file_url = base_url + href
        local_path = os.path.join(download_folder, href)
        if not os.path.exists(local_path):  # Skip if already downloaded
            print(f"Downloading {href}...")
            file_data = requests.get(file_url)
            with open(local_path, 'wb') as f:
                f.write(file_data.content)

print(f"All files downloaded to: {download_folder}")
