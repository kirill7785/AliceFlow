object FormSelProjMat: TFormSelProjMat
  Left = 322
  Top = 114
  Caption = 'Select Project Material'
  ClientHeight = 73
  ClientWidth = 369
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBoxMain: TGroupBox
    Left = 8
    Top = 8
    Width = 361
    Height = 65
    Caption = 'Please, select project material'
    Color = clMoneyGreen
    ParentColor = False
    TabOrder = 0
    object CBlistprojmat: TComboBox
      Left = 8
      Top = 24
      Width = 233
      Height = 21
      TabOrder = 0
    end
    object BApply: TButton
      Left = 264
      Top = 22
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 1
      OnClick = BApplyClick
    end
  end
end
