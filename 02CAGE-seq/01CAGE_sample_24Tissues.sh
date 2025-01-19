#!/bin/bash

###nf-core
workdir=~/FAANG/data/CAGE_data/bigWig
cd $workdir

##modify the download data name to tissues name
#download data from Emily L. Clark group. doi:10.3389/fgene.2020.580580
mv EC1.plus.bw AdrenalCortex.plus.bw
mv EC23.plus.bw AdrenalMedulla.plus.bw
mv EC20.plus.bw Bladder.plus.bw
mv UL33.plus.bw Cerebellum.plus.bw
mv EC3.plus.bw CerebralCortex.plus.bw
mv UL59.plus.bw DescendingColon.plus.bw
mv EC26.plus.bw Duodenum.plus.bw
mv UL31.plus.bw Gallbladder.plus.bw
mv UL62.plus.bw HeartRightAtrium.plus.bw
mv EC45.plus.bw HeartRightVentricle.plus.bw
mv EC7.plus.bw Ileum.plus.bw
mv EC24.plus.bw IleumPeyersPatch.plus.bw
mv UL64.plus.bw Lung.plus.bw
mv EC42.plus.bw LymphNodeMesenteric.plus.bw
mv EC46.plus.bw MuscleSM.plus.bw
mv EC41.plus.bw Ovary.plus.bw
mv EC19.plus.bw Oviduct.plus.bw
mv UL60.plus.bw Reticulum.plus.bw
mv EC37.plus.bw RumenAtrium.plus.bw
mv EC32.plus.bw SpinalCord.plus.bw
mv EC25.plus.bw SpiralColon.plus.bw
mv EC9.plus.bw Tongue.plus.bw
mv EC6.plus.bw Tonsil.plus.bw
mv EC8.plus.bw Uterus.plus.bw

mv EC1.minus.bw AdrenalCortex.minus.bw
mv EC23.minus.bw AdrenalMedulla.minus.bw
mv EC20.minus.bw Bladder.minus.bw
mv UL33.minus.bw Cerebellum.minus.bw
mv EC3.minus.bw CerebralCortex.minus.bw
mv UL59.minus.bw DescendingColon.minus.bw
mv EC26.minus.bw Duodenum.minus.bw
mv UL31.minus.bw Gallbladder.minus.bw
mv UL62.minus.bw HeartRightAtrium.minus.bw
mv EC45.minus.bw HeartRightVentricle.minus.bw
mv EC7.minus.bw Ileum.minus.bw
mv EC24.minus.bw IleumPeyersPatch.minus.bw
mv UL64.minus.bw Lung.minus.bw
mv EC42.minus.bw LymphNodeMesenteric.minus.bw
mv EC46.minus.bw MuscleSM.minus.bw
mv EC41.minus.bw Ovary.minus.bw
mv EC19.minus.bw Oviduct.minus.bw
mv UL60.minus.bw Reticulum.minus.bw
mv EC37.minus.bw RumenAtrium.minus.bw
mv EC32.minus.bw SpinalCord.minus.bw
mv EC25.minus.bw SpiralColon.minus.bw
mv EC9.minus.bw Tongue.minus.bw
mv EC6.minus.bw Tonsil.minus.bw
mv EC8.minus.bw Uterus.minus.bw
