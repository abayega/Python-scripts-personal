awk '//{print($10)}' CCAP.Models_edited_1.gtf2.2 | awk 'BEGIN{FS=";"}//{print($1)}' | sort -u | sed s/"\""/""/g > CCAP.Models_edited_1.gtf2.2_genelist
