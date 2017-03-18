import os
import subprocess


def main():
	sourceImagePath = "../image/"
	targetImagePath = "../equalizedImage/"
	for imageName in os.listdir(sourceImagePath):
		args = ['./Ex3']
		args.append(sourceImagePath + imageName)
		args.append(targetImagePath + imageName)
		subprocess.call(args)


if __name__ == "__main__":
	main()