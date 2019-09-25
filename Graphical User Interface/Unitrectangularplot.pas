unit Unitrectangularplot;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, TeEngine, Series, StdCtrls, ExtCtrls, TeeProcs, Chart;

type
  TfrmRectangularPlot = class(TForm)
    cht1: TChart;
    btnclose: TButton;
    lnsrsSeries1: TLineSeries;
    procedure btncloseClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  frmRectangularPlot: TfrmRectangularPlot;

implementation

{$R *.dfm}

procedure TfrmRectangularPlot.btncloseClick(Sender: TObject);
begin
   Close;
end;

end.
