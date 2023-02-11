object FormAdaptRedionALICEMesh: TFormAdaptRedionALICEMesh
  Left = 0
  Top = 0
  Caption = 'Adapt Region ALICE Mesh'
  ClientHeight = 275
  ClientWidth = 750
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBoxCabinet: TGroupBox
    Left = 8
    Top = 0
    Width = 737
    Height = 49
    Caption = 'Cabinet'
    TabOrder = 0
    object Label1: TLabel
      Left = 16
      Top = 24
      Width = 81
      Height = 13
      Caption = 'Cabinet min level'
    end
    object ComboBoxCabinetMinLevel: TComboBox
      Left = 112
      Top = 21
      Width = 57
      Height = 21
      ItemIndex = 0
      TabOrder = 0
      Text = '-1'
      Items.Strings = (
        '-1'
        '0'
        '1'
        '2'
        '3'
        '4'
        '5'
        '6'
        '7'
        '8'
        '9'
        '10'
        '11'
        '12'
        '13'
        '14'
        '15'
        '16'
        '17'
        '18'
        '19'
        '20'
        '21'
        '22'
        '23'
        '24'
        '25'
        '26'
        '27'
        '28'
        '29'
        '30')
    end
  end
  object GroupBoxAdaptiveRegion1: TGroupBox
    Left = 8
    Top = 55
    Width = 737
    Height = 58
    Caption = 'Adaptive Region 1'
    TabOrder = 1
    object LabelNameminLevel1: TLabel
      Left = 16
      Top = 24
      Width = 41
      Height = 13
      Caption = 'min level'
    end
    object LabelNamexS1: TLabel
      Left = 111
      Top = 24
      Width = 12
      Height = 13
      Caption = 'xS'
    end
    object LabelUnitxS1: TLabel
      Left = 186
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object LabelNamexE1: TLabel
      Left = 219
      Top = 24
      Width = 12
      Height = 13
      Caption = 'xE'
    end
    object LabelUnitxE1: TLabel
      Left = 290
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object LabelNameyS1: TLabel
      Left = 329
      Top = 24
      Width = 12
      Height = 13
      Caption = 'yS'
    end
    object LabelUnityS1: TLabel
      Left = 400
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object LabelNameyE1: TLabel
      Left = 432
      Top = 24
      Width = 12
      Height = 13
      Caption = 'yE'
    end
    object LabelUnityE1: TLabel
      Left = 503
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object LabelNamezS1: TLabel
      Left = 536
      Top = 24
      Width = 11
      Height = 13
      Caption = 'zS'
    end
    object LabelUnitzS1: TLabel
      Left = 612
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object LabelNamezE1: TLabel
      Left = 634
      Top = 24
      Width = 11
      Height = 13
      Caption = 'zE'
    end
    object LabelUnitzE1: TLabel
      Left = 704
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object ComboBoxAdaptiveRegion1: TComboBox
      Left = 63
      Top = 16
      Width = 42
      Height = 21
      ItemIndex = 0
      TabOrder = 0
      Text = '-1'
      Items.Strings = (
        '-1'
        '0'
        '1'
        '2'
        '3'
        '4'
        '5'
        '6'
        '7'
        '8'
        '9'
        '10'
        '11'
        '12'
        '13'
        '14'
        '15'
        '16'
        '17'
        '18'
        '19'
        '20'
        '21'
        '22'
        '23'
        '24'
        '25'
        '26'
        '27'
        '28'
        '29'
        '30')
    end
    object EditxS1: TEdit
      Left = 129
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 1
    end
    object EditxE1: TEdit
      Left = 237
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 2
    end
    object EdityS1: TEdit
      Left = 347
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 3
    end
    object EdityE1: TEdit
      Left = 450
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 4
    end
    object EditzS1: TEdit
      Left = 558
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 5
    end
    object EditzE1: TEdit
      Left = 651
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 6
    end
  end
  object GroupBoxAdaptiveRegion2: TGroupBox
    Left = 5
    Top = 119
    Width = 737
    Height = 58
    Caption = 'Adaptive Region 2'
    TabOrder = 2
    object LabelNameminLevel2: TLabel
      Left = 16
      Top = 24
      Width = 41
      Height = 13
      Caption = 'min level'
    end
    object LabelNamexS2: TLabel
      Left = 111
      Top = 24
      Width = 12
      Height = 13
      Caption = 'xS'
    end
    object LabelUnitxS2: TLabel
      Left = 186
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object LabelNamexE2: TLabel
      Left = 219
      Top = 24
      Width = 12
      Height = 13
      Caption = 'xE'
    end
    object LabelUnitxE2: TLabel
      Left = 290
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object LabelNameyS2: TLabel
      Left = 329
      Top = 24
      Width = 12
      Height = 13
      Caption = 'yS'
    end
    object LabelUnityS2: TLabel
      Left = 400
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object LabelNameyE2: TLabel
      Left = 432
      Top = 24
      Width = 12
      Height = 13
      Caption = 'yE'
    end
    object LabelUnityE2: TLabel
      Left = 503
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object LabelNamezS2: TLabel
      Left = 536
      Top = 24
      Width = 11
      Height = 13
      Caption = 'zS'
    end
    object LabelUnitzS2: TLabel
      Left = 612
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object LabelNamezE2: TLabel
      Left = 634
      Top = 24
      Width = 11
      Height = 13
      Caption = 'zE'
    end
    object LabelUnitzE2: TLabel
      Left = 704
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object ComboBoxAdaptiveRegion2: TComboBox
      Left = 63
      Top = 16
      Width = 42
      Height = 21
      ItemIndex = 0
      TabOrder = 0
      Text = '-1'
      Items.Strings = (
        '-1'
        '0'
        '1'
        '2'
        '3'
        '4'
        '5'
        '6'
        '7'
        '8'
        '9'
        '10'
        '11'
        '12'
        '13'
        '14'
        '15'
        '16'
        '17'
        '18'
        '19'
        '20'
        '21'
        '22'
        '23'
        '24'
        '25'
        '26'
        '27'
        '28'
        '29'
        '30')
    end
    object EditxS2: TEdit
      Left = 129
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 1
    end
    object EditxE2: TEdit
      Left = 237
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 2
    end
    object EdityS2: TEdit
      Left = 347
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 3
    end
    object EdityE2: TEdit
      Left = 450
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 4
    end
    object EditzS2: TEdit
      Left = 558
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 5
    end
    object EditzE2: TEdit
      Left = 651
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 6
    end
  end
  object GroupBoxAdaptiveRegion3: TGroupBox
    Left = 5
    Top = 183
    Width = 737
    Height = 58
    Caption = 'Adaptive Region 3'
    TabOrder = 3
    object LabelNameminLevel3: TLabel
      Left = 16
      Top = 24
      Width = 41
      Height = 13
      Caption = 'min level'
    end
    object LabelNamexS3: TLabel
      Left = 111
      Top = 24
      Width = 12
      Height = 13
      Caption = 'xS'
    end
    object LabelUnitxS3: TLabel
      Left = 186
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object LabelNamexE3: TLabel
      Left = 219
      Top = 24
      Width = 12
      Height = 13
      Caption = 'xE'
    end
    object LabelUnitxE3: TLabel
      Left = 290
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object LabelNameyS3: TLabel
      Left = 329
      Top = 24
      Width = 12
      Height = 13
      Caption = 'yS'
    end
    object LabelUnityS3: TLabel
      Left = 400
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object LabelNameyE3: TLabel
      Left = 432
      Top = 24
      Width = 12
      Height = 13
      Caption = 'yE'
    end
    object LabelUnityE3: TLabel
      Left = 503
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object LabelNamezS3: TLabel
      Left = 536
      Top = 24
      Width = 11
      Height = 13
      Caption = 'zS'
    end
    object LabelUnitzS3: TLabel
      Left = 612
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object LabelNamezE3: TLabel
      Left = 634
      Top = 24
      Width = 11
      Height = 13
      Caption = 'zE'
    end
    object LabelUnitzE3: TLabel
      Left = 704
      Top = 24
      Width = 16
      Height = 13
      Caption = 'mm'
    end
    object ComboBoxAdaptiveRegion3: TComboBox
      Left = 63
      Top = 16
      Width = 42
      Height = 21
      ItemIndex = 0
      TabOrder = 0
      Text = '-1'
      Items.Strings = (
        '-1'
        '0'
        '1'
        '2'
        '3'
        '4'
        '5'
        '6'
        '7'
        '8'
        '9'
        '10'
        '11'
        '12'
        '13'
        '14'
        '15'
        '16'
        '17'
        '18'
        '19'
        '20'
        '21'
        '22'
        '23'
        '24'
        '25'
        '26'
        '27'
        '28'
        '29'
        '30')
    end
    object EditxS3: TEdit
      Left = 129
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 1
    end
    object EditxE3: TEdit
      Left = 237
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 2
    end
    object EdityS3: TEdit
      Left = 347
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 3
    end
    object EdityE3: TEdit
      Left = 450
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 4
    end
    object EditzS3: TEdit
      Left = 558
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 5
    end
    object EditzE3: TEdit
      Left = 651
      Top = 16
      Width = 47
      Height = 21
      TabOrder = 6
    end
  end
  object ButtonApply: TButton
    Left = 624
    Top = 248
    Width = 75
    Height = 25
    Caption = 'Apply'
    TabOrder = 4
    OnClick = ButtonApplyClick
  end
end
