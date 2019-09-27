# DO NOT EDIT
# This makefile makes sure all linkable targets are
# up-to-date with anything they link to
default:
	echo "Do not invoke directly"

# Rules to remove targets that are older than anything to which they
# link.  This forces Xcode to relink the targets from scratch.  It
# does not seem to check these dependencies itself.
PostBuild.Assignment1_bin.Debug:
/Users/daniellu/Documents/Daniel/Grad\ School/NYU/PhD/Classes/Fall\ 2019/Interactive\ Computer\ Graphics/Assignments/CS-GY-6533/Assignment_1/build/Debug/Assignment1_bin:
	/bin/rm -f /Users/daniellu/Documents/Daniel/Grad\ School/NYU/PhD/Classes/Fall\ 2019/Interactive\ Computer\ Graphics/Assignments/CS-GY-6533/Assignment_1/build/Debug/Assignment1_bin


PostBuild.Assignment1_bin.Release:
/Users/daniellu/Documents/Daniel/Grad\ School/NYU/PhD/Classes/Fall\ 2019/Interactive\ Computer\ Graphics/Assignments/CS-GY-6533/Assignment_1/build/Release/Assignment1_bin:
	/bin/rm -f /Users/daniellu/Documents/Daniel/Grad\ School/NYU/PhD/Classes/Fall\ 2019/Interactive\ Computer\ Graphics/Assignments/CS-GY-6533/Assignment_1/build/Release/Assignment1_bin


PostBuild.Assignment1_bin.MinSizeRel:
/Users/daniellu/Documents/Daniel/Grad\ School/NYU/PhD/Classes/Fall\ 2019/Interactive\ Computer\ Graphics/Assignments/CS-GY-6533/Assignment_1/build/MinSizeRel/Assignment1_bin:
	/bin/rm -f /Users/daniellu/Documents/Daniel/Grad\ School/NYU/PhD/Classes/Fall\ 2019/Interactive\ Computer\ Graphics/Assignments/CS-GY-6533/Assignment_1/build/MinSizeRel/Assignment1_bin


PostBuild.Assignment1_bin.RelWithDebInfo:
/Users/daniellu/Documents/Daniel/Grad\ School/NYU/PhD/Classes/Fall\ 2019/Interactive\ Computer\ Graphics/Assignments/CS-GY-6533/Assignment_1/build/RelWithDebInfo/Assignment1_bin:
	/bin/rm -f /Users/daniellu/Documents/Daniel/Grad\ School/NYU/PhD/Classes/Fall\ 2019/Interactive\ Computer\ Graphics/Assignments/CS-GY-6533/Assignment_1/build/RelWithDebInfo/Assignment1_bin




# For each target create a dummy ruleso the target does not have to exist
