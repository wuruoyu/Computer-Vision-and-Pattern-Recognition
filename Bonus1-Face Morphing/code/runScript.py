import os
import subprocess
import imageio

def main():
	sourceImagePath = "../img/1.jpg"
	destImagePath = "../img/2.jpg"
	outImagePath = "../transformProcess/"
	args = ['./Bouns1']
	args.append(sourceImagePath)
	args.append(destImagePath)
	# interval #
	args.append("40")
	#subprocess.call(args)

	imageName = []
	for _imageName in os.listdir(outImagePath):
		imageName.append(_imageName)
	imageName.sort(key=lambda x: int(x.strip(".jpg")))

	images = []
	for _imageName in imageName:
		images.append(imageio.imread(outImagePath + _imageName))
	imageio.mimsave("../animated.gif", images, fps=20)


if __name__ == "__main__":
	main()
