import os
import subprocess


def main():
	sourceImagePath = "../transferSource/Starry night.jpg"
	targetImagePath = "../image/"
	outputImagePath = "../transferImage/"
	for imageName in os.listdir(targetImagePath):
		args = ['./Ex3']
		args.append(sourceImagePath)
		args.append(targetImagePath + imageName)
		args.append(outputImagePath + imageName)
		subprocess.call(args)


if __name__ == "__main__":
	main()