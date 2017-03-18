import os
import subprocess


def main():
	for imageName in os.listdir("../image/"):
		args = ['./Ex3']
		args.append("../image/" + imageName)
		args.append("../equalizedImage/"+imageName)
		subprocess.call(args)


if __name__ == "__main__":
	main()