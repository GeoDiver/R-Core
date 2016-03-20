require_relative 'pool'
require 'json'
require 'fileutils'

class TestOverviewScript
  def initialize(opt)
    @db_path                  = opt[:db_path]
    @rcore_dir                = opt[:rcore_dir]
    @results                  = {}
    @results[:failed]         = []
    @results[:does_not_exist] = []
    @mutex                    = Mutex.new
    @pool                     = Pool.new(opt[:num_threads])
  end

  def start_analysis(top_limit)
    (1..top_limit).each do |i|
      if Dir.exist?(File.join(@db_path, "GDS#{i}"))
        @pool.schedule("GDS#{i}") { |geo_db| run_analysis(geo_db) }
      else
        @results[:does_not_exist] << i
        next
      end
    end
  ensure
    @pool.shutdown
  end

  def print_results
    unless @results[:failed].empty?
      puts 'Failed DBS:'
      puts @results[:failed]
    end

    unless @results[:does_not_exist].empty? && @results[:failed].empty?
      File.open('overview_failures.json', 'w') { |f| f.puts @results.to_json }
    end
  end

  private

  def overview_cmd(db, params)
    "Rscript #{File.expand_path('overview.R', @rcore_dir)}" \
    " --dbrdata #{File.join(@db_path, db, "#{db}.RData")}" \
    " --rundir '#{File.join(@db_path, db, 'overview')}/'" \
    " --analyse 'Boxplot,PCA'" \
    " --accession #{db} --factor \"#{params[:factor]}\"" \
    " --popA \"#{to_comma_delimited_string(params[:groupa])}\"" \
    " --popB \"#{to_comma_delimited_string(params[:groupb])}\"" \
    " --popname1 'Group1' --popname2 'Group2'" \
    ' --dev TRUE'
  end

  def to_comma_delimited_string(arr)
    arr.each { |e| e.gsub!(/(?<!\\),/, '\,') }
    arr.join(',')
  end

  def run_analysis(geo_db)
    params = generate_params(geo_db)
    overview_run_dir = File.join(@db_path, geo_db, 'overview')
    FileUtils.rm_r overview_run_dir if Dir.exist? overview_run_dir
    FileUtils.mkdir overview_run_dir
    system("#{overview_cmd(geo_db, params)}")
    assert_output(geo_db)
  rescue
    assert_output(geo_db)
  end

  def generate_params(db)
    json_file = File.join(@db_path, db, "#{db}.json")
    data = JSON.parse(IO.read(json_file))
    {
      factor: data['Factors'].first[0],
      groupa: data['Factors'].first[1][0..0],
      groupb: data['Factors'].first[1][1..-1]
    }
  end

  def assert_output(db)
    overview_run_dir = File.join(@db_path, db, 'overview')
    boxplot   = File.join(overview_run_dir, 'boxplot.png')
    data_json = File.join(overview_run_dir, 'data.json')
    return if File.exist?(boxplot) && File.exist?(data_json)
    @mutex.synchronize { @results[:failed] << db }
  end
end

opt = {}

opt[:db_path]     = ARGV[0] # Dir with the GDS dir (has to have right structure)
opt[:rcore_dir]   = ARGV[1] # Dir with the RScripts
opt[:num_threads] = ARGV[2].to_i # Number of Threads
top_limit         = ARGV[3].to_i # Number of Analyses to Run

analysis = TestOverviewScript.new(opt)
analysis.start_analysis(top_limit)
analysis.print_results
