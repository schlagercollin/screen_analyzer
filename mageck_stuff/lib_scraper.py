import sys, csv

if __name__ == "__main__":
    with open(sys.argv[1], "r") as original:
        with open(sys.argv[2], "w") as new:
            original_reader = csv.reader(original)
            original_reader.__next__()
            scraped_writer = csv.writer(new)
            for line in original_reader:
                print(line)
                scraped_line = [line[1],line[2],line[7]]
                scraped_writer.writerow(scraped_line)
