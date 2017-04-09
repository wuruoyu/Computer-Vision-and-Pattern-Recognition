import os
import subprocess


def main():
	sourceImagePath = "../transferSource/"
	targetImagePath = "../image/"
	outputImagePath = "../transferImage/"
	for targetImageName in os.listdir(targetImagePath):
		for sourceImageName in os.listdir(sourceImagePath):
			args = ['./Ex3']
			args.append(sourceImagePath + sourceImageName)
			args.append(targetImagePath + targetImageName)
			args.append(outputImagePath + targetImageName.strip('.jpg') + "_" + sourceImageName)
			subprocess.call(args)


if __name__ == "__main__":
	main()

