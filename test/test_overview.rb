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
    @threads                  = opt[:num_threads]
  end

  def start_analysis(top_limit, type)
    (1..top_limit).each do |i|
      if Dir.exist?(File.join(@db_path, "GDS#{i}"))
        @pool.schedule("GDS#{i}") { |geo_db| run_analysis(geo_db, type) }
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

  def run_analysis(geo_db, type)
    params = generate_params(geo_db)
    run_dir = File.join(@db_path, geo_db, 'run_dir')
    FileUtils.rm_r run_dir if Dir.exist? run_dir
    FileUtils.mkdir run_dir
    if @threads == 1
      STDERR.puts
      STDERR.puts '############'
      STDERR.puts '############'
      STDERR.puts
      STDERR.puts overview_cmd(geo_db, params) if type == 'overview' || type == 'all'
      STDERR.puts dgea_cmd(geo_db, params) if type == 'dgea' || type == 'all'
      STDERR.puts overview_cmd(geo_db, params) if type == 'gage' || type == 'all'
      STDERR.puts
    end
    system("#{overview_cmd(geo_db, params)}") if type == 'overview' || type == 'all'
    system("#{dgea_cmd(geo_db, params)}") if type == 'dgea' || type == 'all'
    system("#{gage_cmd(geo_db, params)}") if type == 'gage' || type == 'all'
  ensure
    assert_overview_output(geo_db) if type == 'overview' || type == 'all'
    assert_dgea_output(geo_db) if type == 'dgea' || type == 'all'
    assert_gage_output(geo_db) if type == 'gage' || type == 'all'
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

  def overview_cmd(db, params)
    "Rscript #{File.expand_path('overview.R', @rcore_dir)}" \
    " --dbrdata #{File.join(@db_path, db, "#{db}.RData")}" \
    " --rundir '#{File.join(@db_path, db, 'run_dir')}/'" \
    " --analyse 'Boxplot,PCA'" \
    " --accession #{db} --factor \"#{params[:factor]}\"" \
    " --popA \"#{to_comma_delimited_string(params[:groupa])}\"" \
    " --popB \"#{to_comma_delimited_string(params[:groupb])}\"" \
    " --popname1 'Group1' --popname2 'Group2'" \
    ' --dev TRUE'
  end

  def dgea_cmd(db, params)
    "Rscript #{File.expand_path('dgea.R', @rcore_dir)}" \
    " --dbrdata #{File.join(@db_path, db, "#{db}.RData")}" \
    " --rundir '#{File.join(@db_path, db, 'run_dir')}/'" \
    " --analyse 'Boxplot,PCA,Volcano,Heatmap'" \
    " --accession #{db} --factor \"#{params[:factor]}\"" \
    " --popA \"#{to_comma_delimited_string(params[:groupa])}\"" \
    " --popB \"#{to_comma_delimited_string(params[:groupb])}\"" \
    " --popname1 'Group1' --popname2 'Group2' --topgenecount 250" \
    " --foldchange 0 --thresholdvalue 0 --distance 'euclidean'" \
    " --clustering 'average' --clusterby 'Complete' --heatmaprows 100 " \
    " --adjmethod 'fdr' --dendrow TRUE --dendcol TRUE --dev TRUE"
  end


  def gage_cmd(db, params)
    "Rscript  #{File.expand_path('gage.R', @rcore_dir)}" \
    " --dbrdata #{File.join(@db_path, db, "#{db}.RData")}" \
    " --rundir '#{File.join(@db_path, db, 'run_dir')}/'" \
    " --accession #{db} --factor \"#{params[:factor]}\"" \
    " --popA \"#{to_comma_delimited_string(params[:groupa])}\"" \
    " --popB \"#{to_comma_delimited_string(params[:groupb])}\"" \
    " --comparisontype 'ExpVsCtrl' --genesettype 'KEGG'" \
    " --distance 'euclidean' --clustering 'average' --clusterby 'Complete'" \
    " --heatmaprows 100 --dendrow TRUE --dendcol TRUE  --dev TRUE"
  end

  def to_comma_delimited_string(arr)
    arr.each do |e|
      e.gsub!(/(?<!\\),/, '\,')
      e.gsub!('-', '\-')
    end
    arr.join(',')
  end

  def assert_overview_output(db)
    run_dir = File.join(@db_path, db, 'run_dir')
    boxplot   = File.join(run_dir, 'boxplot.png')
    data_json = File.join(run_dir, 'data.json')
    return if File.exist?(boxplot) && File.exist?(data_json)
    @mutex.synchronize { @results[:failed] << db }
  end

  def assert_dgea_output(db)
    run_dir      = File.join(@db_path, db, 'run_dir')
    heatmap      = File.join(run_dir, 'dgea_heatmap.svg')
    volcano      = File.join(run_dir, 'dgea_volcano.png')
    toptable     = File.join(run_dir, 'dgea_toptable.RData')
    toptable_tsv = File.join(run_dir, 'dgea_toptable.tsv')
    data_json    = File.join(run_dir, 'dgea_data.json')
    return if File.exist?(heatmap) && File.exist?(volcano) &&
    File.exist?(toptable) && File.exist?(toptable_tsv) && File.exist?(data_json)
    @mutex.synchronize { @results[:failed] << db }
  end

  def assert_gage_output(db)
    run_dir      = File.join(@db_path, db, 'run_dir')
    heatmap      = File.join(run_dir, 'gage_heatmap.svg')
    toptable     = File.join(run_dir, 'gage.RData')
    data_json    = File.join(run_dir, 'gage_data.json')
    toptable_tsv = File.join(run_dir, 'gage_toptable.tsv')
    return if File.exist?(heatmap) && File.exist?(volcano) &&
    File.exist?(toptable) && File.exist?(toptable_tsv) && File.exist?(data_json)
    @mutex.synchronize { @results[:failed] << db }
  end
end

opt = {}

opt[:db_path]     = ARGV[0] # Dir with the GDS dir (has to have right structure)
opt[:rcore_dir]   = ARGV[1] # Dir with the RScripts
opt[:num_threads] = ARGV[2].to_i # Number of Threads
top_limit         = ARGV[3].to_i # Number of Analyses to Run

type = 'all' # overview, dgea, gage, all 

analysis = TestOverviewScript.new(opt)
analysis.start_analysis(top_limit, type)
analysis.print_results
