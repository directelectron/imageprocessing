#!/usr/bin/env python3
import os
import sys

def help():
	print("Usage: python mdoc_xml.py <command>")
	print("Modified from https://github.com/ccgauvin94/cs_jiffies/blob/main/mdoc_xml.py")
	print("\nCommands:")
	print(" --help, -h	  -- Display help")
	print(" --dir,  -d	  -- Specify directory")
	print(" --mrc,  -m	  -- Use .mrc.xml for output instead of .tif.xml")
	print("\nWill default to current directory unless specified.")

def get_mdocs(directory):
	mdocs = []
	with os.scandir(directory) as entries:
		for entry in entries:
			if entry.is_file() and entry.name.endswith(".mdoc"):
				mdocs.append(entry.name)
	print("Found", len(mdocs), "mdoc files in", str(directory))
	return mdocs

def extract_image_shift(mdoc):
	image_shift_x = ""
	image_shift_y = ""
	with open(mdoc, 'r') as file:
		for line in file:
			if "ImageShift" in line:
				image_shift = line.strip()
				image_shift_x = (str(image_shift.rsplit()[2]))
				image_shift_y = (str(image_shift.rsplit()[3]))
				return(image_shift_x, image_shift_y)
				break

def write_xml(mdoc, image_ext):
	is_x, is_y = extract_image_shift(mdoc)
	print(str(mdoc), "has image shift of", is_x, is_y)
	xml_name = mdoc.rsplit(".mdoc",2)[0].rsplit(".tif",2)[0].rsplit(".mrc",2)[0]
	xml_file = open((str(xml_name) + image_ext + ".xml"), "w")
	xml_contents = """<MicroscopeImage xmlns="http://schemas.datacontract.org/2004/07/Fei.SharedObjects" xmlns:i="http://www.w3.org/2001/XMLSchema-instance">
	<microscopeData>
		<optics>
			<BeamShift xmlns:a="http://schemas.datacontract.org/2004/07/Fei.Types">
				<a:_x>%s</a:_x>
				<a:_y>%s</a:_y>
			</BeamShift>
		</optics>
	</microscopeData>
</MicroscopeImage>
""" % (is_x, is_y)
	xml_file.writelines(xml_contents)
	xml_file.close()
	return

def main():
	
	directory = "."
	image_ext = ".tif"
	
	if len(sys.argv) > 1:
		i = 1
		while i < len(sys.argv):
			if (sys.argv[i] == "--help") or (sys.argv[i] == "-h"):
				help()
				sys.exit(0)
			elif (sys.argv[i] == "--dir") or (sys.argv[i] == "-d"):
				i += 1
				directory = sys.argv[i]
			elif (sys.argv[i] == "--mrc") or (sys.argv[i] == "-m"):
				image_ext = ".mrc"
			i += 1
	
	mdoc_files = get_mdocs(directory)
	for i in range (len(mdoc_files)):
		write_xml(mdoc_files[i], image_ext)

if __name__ == "__main__":
	main()
