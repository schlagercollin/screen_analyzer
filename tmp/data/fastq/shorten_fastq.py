def shorten_fastq(file_name, shortened_name, sequences):
    with open(file_name, 'r') as long_file:
        with open(shortened_name, 'w') as short_file:
            for count, line in enumerate(long_file):
                if count >= (sequences*4):
                    break
                else:
                    short_file.write(line)
