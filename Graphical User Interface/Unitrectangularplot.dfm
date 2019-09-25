object frmRectangularPlot: TfrmRectangularPlot
  Left = 373
  Top = 193
  Width = 779
  Height = 602
  Caption = 'Rectangular Plot'
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object cht1: TChart
    Left = 0
    Top = 0
    Width = 761
    Height = 529
    BackWall.Brush.Color = clWhite
    BackWall.Brush.Style = bsClear
    Title.Text.Strings = (
      'TChart')
    View3D = False
    View3DWalls = False
    TabOrder = 0
    object lnsrsSeries1: TLineSeries
      Marks.ArrowLength = 8
      Marks.Visible = False
      SeriesColor = clRed
      Pointer.InflateMargins = True
      Pointer.Style = psRectangle
      Pointer.Visible = False
      XValues.DateTime = False
      XValues.Name = 'X'
      XValues.Multiplier = 1.000000000000000000
      XValues.Order = loAscending
      YValues.DateTime = False
      YValues.Name = 'Y'
      YValues.Multiplier = 1.000000000000000000
      YValues.Order = loNone
    end
  end
  object btnclose: TButton
    Left = 656
    Top = 536
    Width = 75
    Height = 25
    Caption = 'Close'
    TabOrder = 1
    OnClick = btncloseClick
  end
end
