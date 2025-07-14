import os
import re

# Putanja do foldera sa fajlovima
folder_path = "asteroids"

# RegEx za pronalaženje decimalnih brojeva u imenu fajla
pattern = re.compile(r"(.*?_)(\d+\.\d+)([_].*)")

for filename in os.listdir(folder_path):
    old_path = os.path.join(folder_path, filename)

    # Preskoči ako nije fajl
    if not os.path.isfile(old_path):
        continue

    match = pattern.match(filename)
    if match:
        prefix, number_str, suffix = match.groups()
        try:
            number = round(float(number_str), 1)
            # Formatiraj broj tako da uvek ima jednu decimalu (npr. 2.0)
            number_formatted = f"{number:.1f}"
            new_filename = f"{prefix}{number_formatted}{suffix}"
            new_path = os.path.join(folder_path, new_filename)

            if new_filename != filename:
                os.rename(old_path, new_path)
                print(f"Renamed: {filename} -> {new_filename}")
        except ValueError:
            continue
