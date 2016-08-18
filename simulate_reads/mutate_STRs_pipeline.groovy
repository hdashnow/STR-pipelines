
mutate_background = {
    produce('background.bed') {
        exec """
            /Users/hd_vlsci/Documents/git/microsat_stats/STR_simulation_script.R
                -L chr2:233712201-233712246
                /Users/hd_vlsci/Documents/reference-data/hg19.simpleRepeat.txt.gz
                /Users/hd_vlsci/Documents/reference-data/str-stats
                --background $output.bed
        """

    }
}

mutate_locus = {
    branch.simID = branch.name

    produce(branch.simID + '.bed') {
        exec """
            /Users/hd_vlsci/Documents/git/microsat_stats/STR_simulation_script.R
                -L chr2:233712201-233712246
                /Users/hd_vlsci/Documents/reference-data/hg19.simpleRepeat.txt.gz
                /Users/hd_vlsci/Documents/reference-data/str-stats
                > $output.bed
        """

    }
}

simID = (1..10)

run {
    mutate_background +
    simID * [mutate_locus]

}
